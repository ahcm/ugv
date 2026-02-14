use anyhow::Result;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

/// A single data point in a TSV track
#[derive(Debug, Clone)]
pub struct TsvPoint
{
    pub chromosome: String,
    pub start: usize,
    pub end: usize,
    pub label: String,
}

/// Track for one chromosome
#[derive(Clone)]
pub struct TsvChromosomeTrack
{
    pub points: Vec<TsvPoint>,
}

impl TsvChromosomeTrack
{
    pub fn new() -> Self
    {
        Self { points: Vec::new() }
    }

    /// Add a point and update min/max
    pub fn add_point(&mut self, point: TsvPoint)
    {
        self.points.push(point);
    }

    /// Sort points by position
    pub fn sort(&mut self)
    {
        self.points.sort_by_key(|p| p.start);
    }

    /// Get points in a region
    pub fn get_points_in_region(&self, start: usize, end: usize) -> Vec<&TsvPoint>
    {
        self.points
            .iter()
            .filter(|p| p.end >= start && p.start <= end)
            .collect()
    }
}

/// Container for all TSV tracks
#[derive(Clone)]
pub struct TsvData
{
    pub tracks: HashMap<String, TsvChromosomeTrack>,
    pub track_name: String,
}

impl TsvData
{
    /// Create new TSV data container
    pub fn new(track_name: String) -> Self
    {
        Self {
            tracks: HashMap::new(),
            track_name,
        }
    }

    /// Parse TSV from file (native only)
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_file(path: &str) -> Result<Self>
    {
        use std::fs::File;
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        let track_name = std::path::Path::new(path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("Custom Track")
            .to_string();

        Self::parse_tsv(reader, track_name)
    }

    /// Parse TSV from bytes
    pub fn from_bytes(data: Vec<u8>, track_name: String) -> Result<Self>
    {
        let reader = BufReader::new(&data[..]);
        Self::parse_tsv(reader, track_name)
    }

    /// Parse TSV data
    fn parse_tsv<R: BufRead>(reader: R, track_name: String) -> Result<Self>
    {
        let mut tsv_data = Self::new(track_name);

        for (line_num, line_result) in reader.lines().enumerate()
        {
            let line = line_result?;
            let line = line.trim();

            // Skip empty lines and comments
            if line.is_empty() || line.starts_with('#')
            {
                continue;
            }

            // Parse tab-separated fields
            let fields: Vec<&str> = line.split('\t').collect();

            if fields.len() < 4
            {
                eprintln!(
                    "Warning: Line {} has fewer than 4 fields, skipping: {}",
                    line_num + 1,
                    line
                );
                continue;
            }

            let chromosome = fields[0].to_string();

            let start = match fields[1].parse::<usize>()
            {
                Ok(pos) => pos,
                Err(e) =>
                {
                    eprintln!(
                        "Warning: Line {} has invalid start position '{}': {}",
                        line_num + 1,
                        fields[1],
                        e
                    );
                    continue;
                }
            };

            let end = match fields[2].parse::<usize>()
            {
                Ok(pos) => pos,
                Err(e) =>
                {
                    eprintln!(
                        "Warning: Line {} has invalid end position '{}': {}",
                        line_num + 1,
                        fields[2],
                        e
                    );
                    continue;
                }
            };

            let label = fields[3].to_string();

            // Add point to appropriate chromosome track
            let track = tsv_data
                .tracks
                .entry(chromosome.clone())
                .or_insert_with(TsvChromosomeTrack::new);

            track.add_point(TsvPoint {
                chromosome,
                start,
                end,
                label,
            });
        }

        // Sort all tracks by position
        for track in tsv_data.tracks.values_mut()
        {
            track.sort();
        }

        Ok(tsv_data)
    }
}

#[cfg(test)]
mod tests
{
    use super::*;

    #[test]
    fn test_parse_tsv()
    {
        let data =
            b"chr1\t100\t150\tfeature1\nchr1\t200\t250\tfeature2\nchr2\t150\t180\tfeature3\n";
        let tsv = TsvData::from_bytes(data.to_vec(), "test".to_string()).unwrap();

        assert_eq!(tsv.tracks.len(), 2);
        assert!(tsv.tracks.contains_key("chr1"));
        assert!(tsv.tracks.contains_key("chr2"));

        let chr1_track = &tsv.tracks["chr1"];
        assert_eq!(chr1_track.points.len(), 2);
        assert_eq!(chr1_track.points[0].start, 100);
        assert_eq!(chr1_track.points[0].end, 150);
        assert_eq!(chr1_track.points[0].label, "feature1");
        assert_eq!(chr1_track.points[1].start, 200);
        assert_eq!(chr1_track.points[1].end, 250);
        assert_eq!(chr1_track.points[1].label, "feature2");
    }

    #[test]
    fn test_get_points_in_region()
    {
        let mut track = TsvChromosomeTrack::new();
        track.add_point(TsvPoint {
            chromosome: "chr1".to_string(),
            start: 100,
            end: 150,
            label: "feature1".to_string(),
        });
        track.add_point(TsvPoint {
            chromosome: "chr1".to_string(),
            start: 200,
            end: 250,
            label: "feature2".to_string(),
        });
        track.add_point(TsvPoint {
            chromosome: "chr1".to_string(),
            start: 300,
            end: 320,
            label: "feature3".to_string(),
        });

        let points = track.get_points_in_region(150, 250);
        assert_eq!(points.len(), 2);
        assert_eq!(points[0].start, 100);
        assert_eq!(points[1].start, 200);
    }
}
