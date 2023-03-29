/*!

# PlasMod - apply modula function to alignments to account for circular nature of plasmid

## Usage

Before alignment, the plasmid reference should be duplicated to account for reads which bridge
the "end" of the plasmid.  A standard aligner (bwa, minimap2) is then used to align the reads
to the duplicated reference, after which plasmod translates the BAM file into unduplicated
coordinates.


Eg, a plasmid with sequence AATTGGCC and read GGCCAATT will yield this alignment (space inserted
for emphasis):

AATTGGCCAATTGGCC
----GGCCaatt----

Plasmod will then translate that alignment into the original reference coordinates:

AATTGGCC
aattGGCC

Gaps are represented by RefSkip, unless
Note that only primary alignments are preserved; output is written to stdout.
 */

extern crate log;
use anyhow::Result;
use log::debug;
use rust_htslib::bam::record::{Cigar, CigarString, Record};
use rust_htslib::bam::{Format, Header, Read, Reader, Writer};
use rust_htslib::htslib;
use std::error::Error;
use std::iter;
use std::path::PathBuf;

pub const NONPRIMARY: u16 = (htslib::BAM_FUNMAP
    | htslib::BAM_FSECONDARY
    | htslib::BAM_FQCFAIL
    | htslib::BAM_FDUP
    | htslib::BAM_FSUPPLEMENTARY) as u16;

pub fn halve(ref_len: usize, bam_path: &PathBuf, use_del: bool) -> Result<()> {
    let reader: Reader = Reader::from_path(bam_path).unwrap();
    let hdr = reader.header();

    let mut writer = Writer::from_stdout(&Header::from_template(&hdr), Format::Sam)?;

    for aln in Reader::from_path(bam_path)
        .unwrap()
        .rc_records()
        .map(|x| x.unwrap())
    {
        if aln.flags() & NONPRIMARY == 0 {
            let (new_start, mapped_cigar) =
                mod_cigar(ref_len, aln.pos() as usize, aln.cigar().iter(), use_del);

            let mut new_aln = Record::new();
            new_aln.set(
                aln.qname(),
                Some(&CigarString(mapped_cigar)),
                aln.seq().as_bytes().as_slice(),
                aln.qual(),
            );
            new_aln.set_pos(new_start.try_into().unwrap());

            writer.write(&new_aln)?;
        }
    }

    Ok(())
}

/// return a new cigar w/ the same type as ref_cig
fn new_cigar(ref_cig: &Cigar, len: u32) -> Cigar {
    match ref_cig {
        Cigar::Match(_) => Cigar::Match(len),
        Cigar::Ins(_) => Cigar::Ins(len),
        Cigar::Del(_) => Cigar::Del(len),
        Cigar::RefSkip(_) => Cigar::RefSkip(len),
        Cigar::SoftClip(_) => Cigar::SoftClip(len),
        Cigar::HardClip(_) => Cigar::HardClip(len),
        Cigar::Pad(_) => Cigar::Pad(len),
        Cigar::Equal(_) => Cigar::Equal(len),
        Cigar::Diff(_) => Cigar::Diff(len),
    }
}

/// length of reference occupied by cigar string - eg, Ins(N) doesn't consume any reference
fn ref_occupancy<'a>(mut cigars: std::slice::Iter<'a, Cigar>) -> u32 {
    cigars
        .map(|c| match c {
            Cigar::Match(len) => *len,
            Cigar::Del(len) => *len,
            Cigar::RefSkip(len) => *len,
            Cigar::Equal(len) => *len,
            Cigar::Diff(len) => *len,
            _ => 0,
        })
        .sum()
}

/// apply "modula" to cigar string
fn mod_cigar<'a>(
    ref_len: usize,
    aln_pos: usize,
    mut cigars: std::slice::Iter<'a, Cigar>,
    use_del: bool,
) -> (usize, Vec<Cigar>) {
    let mut pos: usize = aln_pos % ref_len;

    // if the alignment starts in the second copy, then the cigar string can
    // be used w/out modification, and only the start position needs to change
    if aln_pos >= ref_len {
        return (pos, cigars.cloned().collect());
    }

    // for operations that land in the first reference copy
    let mut suffix: Vec<Cigar> = vec![];

    // for operations that wrap around into the second reference copy
    let mut prefix: Vec<Cigar> = vec![];

    // we want to divide the cigar operations, stashing those from the first half
    // in @suffix, and those from the second half in @prefix
    for op in cigars {
        if pos >= ref_len {
            prefix.push(op.clone());
        } else if (pos + op.len() as usize) < ref_len {
            suffix.push(op.clone());
        } else {
            // this operation spans the boundary; cut it
            let first_half = ref_len - pos;
            let second_half = op.len() as usize - first_half;
            suffix.push(new_cigar(&op, first_half as u32));
            prefix.push(new_cigar(&op, second_half as u32));
        }
        pos += ref_occupancy(vec![op.clone()].iter()) as usize;
    }

    if prefix.len() == 0 {
        // trivial case - everything landed in the first copy of the reference
        return (aln_pos % ref_len, suffix);
    }

    let end_of_prefix = ref_occupancy(prefix.iter());
    let len_of_suffix = ref_occupancy(suffix.iter());
    let gap_len =
        <usize as TryInto<u32>>::try_into(ref_len).unwrap() - end_of_prefix - len_of_suffix;
    let gap = if use_del {
        Cigar::Del(gap_len)
    } else {
        Cigar::RefSkip(gap_len)
    };
    (
        0,
        prefix
            .into_iter()
            .chain(iter::once(gap))
            .chain(suffix.into_iter())
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simplest_case() {
        // ref seq is 100bp; alignment is simple 10bp match at start
        assert_eq!(
            (0, vec![Cigar::Match(10)]),
            mod_cigar(100, 0, vec![Cigar::Match(10)].iter(), true)
        );
    }

    #[test]
    fn test_entirely_second() {
        // ref seq is 100bp; alignment is simple 10bp match starting at position 100,
        // ie, entirely in the second duplicate
        assert_eq!(
            (0, vec![Cigar::Match(10)]),
            mod_cigar(100, 100, vec![Cigar::Match(10)].iter(), true)
        );
    }

    #[test]
    fn test_spanning() {
        assert_eq!(
            (
                0,
                vec![Cigar::Match(10), Cigar::RefSkip(80), Cigar::Match(10)]
            ),
            mod_cigar(100, 90, vec![Cigar::Match(20)].iter(), false)
        );
    }

    #[test]
    fn test_spanning_del() {
        assert_eq!(
            (0, vec![Cigar::Match(10), Cigar::Del(80), Cigar::Match(10)]),
            mod_cigar(100, 90, vec![Cigar::Match(20)].iter(), true)
        );
    }

    #[test]
    fn test_mix() {
        assert_eq!(
            (
                0,
                vec![
                    Cigar::Match(10),
                    Cigar::RefSkip(40),
                    Cigar::Match(20),
                    Cigar::Ins(20),
                    Cigar::Match(30)
                ]
            ),
            mod_cigar(
                100,
                50,
                vec![Cigar::Match(20), Cigar::Ins(20), Cigar::Match(40)].iter(),
                false
            )
        );
    }
}
