/// Codon to amino acid translation table (standard genetic code)
pub fn translate_codon(codon: &[u8]) -> char
{
    if codon.len() != 3
    {
        return 'X';
    }

    match codon
    {
        // TTT, TTC -> F (Phe)
        b"TTT" | b"TTC" => 'F',
        // TTA, TTG, CTT, CTC, CTA, CTG -> L (Leu)
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        // ATT, ATC, ATA -> I (Ile)
        b"ATT" | b"ATC" | b"ATA" => 'I',
        // ATG -> M (Met)
        b"ATG" => 'M',
        // GTT, GTC, GTA, GTG -> V (Val)
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        // TCT, TCC, TCA, TCG, AGT, AGC -> S (Ser)
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
        // CCT, CCC, CCA, CCG -> P (Pro)
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        // ACT, ACC, ACA, ACG -> T (Thr)
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        // GCT, GCC, GCA, GCG -> A (Ala)
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        // TAT, TAC -> Y (Tyr)
        b"TAT" | b"TAC" => 'Y',
        // TAA, TAG, TGA -> * (Stop)
        b"TAA" | b"TAG" | b"TGA" => '*',
        // CAT, CAC -> H (His)
        b"CAT" | b"CAC" => 'H',
        // CAA, CAG -> Q (Gln)
        b"CAA" | b"CAG" => 'Q',
        // AAT, AAC -> N (Asn)
        b"AAT" | b"AAC" => 'N',
        // AAA, AAG -> K (Lys)
        b"AAA" | b"AAG" => 'K',
        // GAT, GAC -> D (Asp)
        b"GAT" | b"GAC" => 'D',
        // GAA, GAG -> E (Glu)
        b"GAA" | b"GAG" => 'E',
        // TGT, TGC -> C (Cys)
        b"TGT" | b"TGC" => 'C',
        // TGG -> W (Trp)
        b"TGG" => 'W',
        // CGT, CGC, CGA, CGG, AGA, AGG -> R (Arg)
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
        // GGT, GGC, GGA, GGG -> G (Gly)
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
        // Unknown or ambiguous
        _ => 'X',
    }
}

/// Translate DNA sequence to amino acids in a specific reading frame
/// frame: 0, 1, or 2 for the offset
pub fn translate_frame(sequence: &[u8], frame: usize) -> Vec<char>
{
    let mut amino_acids = Vec::new();

    let mut i = frame;
    while i + 2 < sequence.len()
    {
        let codon = &sequence[i..i + 3];
        amino_acids.push(translate_codon(codon));
        i += 3;
    }

    amino_acids
}

/// Get reverse complement of a DNA sequence
pub fn reverse_complement(sequence: &[u8]) -> Vec<u8>
{
    sequence
        .iter()
        .rev()
        .map(|&base| complement_base(base))
        .collect()
}

/// Get complement of a single base
fn complement_base(base: u8) -> u8
{
    match base
    {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'a' => b't',
        b't' => b'a',
        b'g' => b'c',
        b'c' => b'g',
        _ => b'N',
    }
}

/// Translate all 6 reading frames
/// Returns: (forward frames [0,1,2], reverse frames [0,1,2])
pub fn translate_six_frames(sequence: &[u8]) -> (Vec<Vec<char>>, Vec<Vec<char>>)
{
    let forward_frames = vec![
        translate_frame(sequence, 0),
        translate_frame(sequence, 1),
        translate_frame(sequence, 2),
    ];

    let rev_comp = reverse_complement(sequence);
    let reverse_frames = vec![
        translate_frame(&rev_comp, 0),
        translate_frame(&rev_comp, 1),
        translate_frame(&rev_comp, 2),
    ];

    (forward_frames, reverse_frames)
}

#[cfg(test)]
mod tests
{
    use super::*;

    #[test]
    fn test_translate_codon()
    {
        assert_eq!(translate_codon(b"ATG"), 'M'); // Start codon
        assert_eq!(translate_codon(b"TAA"), '*'); // Stop codon
        assert_eq!(translate_codon(b"GGG"), 'G'); // Glycine
        assert_eq!(translate_codon(b"TTT"), 'F'); // Phenylalanine
    }

    #[test]
    fn test_reverse_complement()
    {
        let seq = b"ATCG";
        let rev_comp = reverse_complement(seq);
        assert_eq!(rev_comp, b"CGAT");
    }

    #[test]
    fn test_translate_frame()
    {
        let seq = b"ATGGGGAAATAA"; // ATG GGG AAA TAA -> M G K *
        let amino_acids = translate_frame(seq, 0);
        assert_eq!(amino_acids, vec!['M', 'G', 'K', '*']);
    }
}
