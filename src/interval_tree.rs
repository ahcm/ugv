use crate::gff::Feature;
use std::collections::HashMap;

#[derive(Debug)]
pub struct IntervalTree
{
    trees: HashMap<String, ChromosomeTree>,
}

#[derive(Debug)]
struct ChromosomeTree
{
    intervals: Vec<Interval>,
}

#[derive(Debug, Clone)]
struct Interval
{
    start: usize,
    end: usize,
    feature_idx: usize,
}

impl IntervalTree
{
    pub fn from_features(features: &[Feature]) -> Self
    {
        let mut trees: HashMap<String, Vec<Interval>> = HashMap::new();

        for (idx, feature) in features.iter().enumerate()
        {
            trees
                .entry(feature.seqid.clone())
                .or_default()
                .push(Interval {
                    start: feature.start,
                    end: feature.end,
                    feature_idx: idx,
                });
        }

        let mut interval_tree = IntervalTree {
            trees: HashMap::new(),
        };

        for (chr, mut intervals) in trees
        {
            intervals.sort_by_key(|i| i.start);
            interval_tree
                .trees
                .insert(chr, ChromosomeTree { intervals });
        }

        interval_tree
    }

    pub fn query(&self, chr: &str, start: usize, end: usize) -> Vec<usize>
    {
        if let Some(tree) = self.trees.get(chr)
        {
            tree.query(start, end)
        }
        else
        {
            Vec::new()
        }
    }
}

impl ChromosomeTree
{
    fn query(&self, start: usize, end: usize) -> Vec<usize>
    {
        let mut result = Vec::new();

        // Binary search to find the first potentially overlapping interval
        let mut left = 0;
        let mut right = self.intervals.len();

        while left < right
        {
            let mid = (left + right) / 2;
            if self.intervals[mid].end <= start
            {
                left = mid + 1;
            }
            else
            {
                right = mid;
            }
        }

        // Collect all overlapping intervals from this point
        for interval in &self.intervals[left..]
        {
            if interval.start >= end
            {
                break;
            }
            if interval.end > start
            {
                result.push(interval.feature_idx);
            }
        }

        result
    }
}

#[cfg(test)]
mod tests
{
    use super::*;

    #[test]
    fn test_interval_query()
    {
        let intervals = vec![
            Interval {
                start: 100,
                end: 200,
                feature_idx: 0,
            },
            Interval {
                start: 150,
                end: 250,
                feature_idx: 1,
            },
            Interval {
                start: 300,
                end: 400,
                feature_idx: 2,
            },
        ];

        let tree = ChromosomeTree { intervals };

        let result = tree.query(120, 180);
        assert_eq!(result, vec![0, 1]);

        let result = tree.query(250, 350);
        assert_eq!(result, vec![2]);

        let result = tree.query(500, 600);
        assert_eq!(result.len(), 0);
    }
}
