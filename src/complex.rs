use std::collections::BTreeSet;

#[derive(Debug)]
pub struct Complex {
    pub(crate) proteins: BTreeSet<String>,
}

impl Complex {
    // 计算两个复合物的匹配系数和重叠的蛋白质个数
    pub fn os(&self, other: &Self) -> (f64, u32) {
        if self.is_empty() || other.is_empty() {
            return (0.0, 0);
        }

        let common = self.proteins.intersection(&other.proteins).count() as f64;
        let os_score = common.powf(2.0) / ((self.size() * other.size()) as f64);
        (os_score, common as u32)
    }

    pub fn size(&self) -> u16 {
        return self.proteins.len() as u16;
    }

    pub fn is_empty(&self) -> bool {
        return self.proteins.is_empty();
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;
    use super::Complex;

    #[test]
    fn test_os() {
        let a = Complex {
            proteins: BTreeSet::from(["A".to_string(), "B".to_string(), "C".to_string()]),
        };
        let b = Complex {
            proteins: BTreeSet::from(["A".to_string(), "B".to_string(), "D".to_string()]),
        };
        let os_ab = a.os(&b);
        println!("{:?}", os_ab);
    }
}
