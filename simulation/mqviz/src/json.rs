//! Minimal JSON value, printer and parser (no dependencies): run records and saved states.

#[derive(Clone, Debug, PartialEq)]
pub enum J { Null, Bool(bool), Num(f64), Str(String), Arr(Vec<J>), Obj(Vec<(String, J)>) }

impl J {
    pub fn obj() -> J { J::Obj(vec![]) }
    /// Builder: add a key.
    pub fn with(mut self, k: &str, v: impl Into<J>) -> J { if let J::Obj(o) = &mut self { o.push((k.to_string(), v.into())); } self }
    pub fn set(&mut self, k: &str, v: impl Into<J>) { if let J::Obj(o) = self { o.push((k.to_string(), v.into())); } }
    pub fn get(&self, k: &str) -> Option<&J> { if let J::Obj(o) = self { o.iter().find(|(kk, _)| kk == k).map(|(_, v)| v) } else { None } }
    pub fn f(&self, k: &str) -> Option<f64> { match self.get(k)? { J::Num(x) => Some(*x), _ => None } }
    pub fn s(&self, k: &str) -> Option<&str> { match self.get(k)? { J::Str(x) => Some(x), _ => None } }
    pub fn b(&self, k: &str) -> Option<bool> { match self.get(k)? { J::Bool(x) => Some(*x), _ => None } }
    pub fn arr(&self, k: &str) -> Option<&Vec<J>> { match self.get(k)? { J::Arr(x) => Some(x), _ => None } }

    pub fn pretty(&self) -> String { let mut s = String::new(); self.write(&mut s, 0); s.push('\n'); s }

    fn write(&self, out: &mut String, ind: usize) {
        let pad = |n: usize| "  ".repeat(n);
        match self {
            J::Null => out.push_str("null"),
            J::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            J::Num(x) => {
                if !x.is_finite() { out.push_str("null"); }
                else if x.fract() == 0.0 && x.abs() < 1e15 { out.push_str(&format!("{}", *x as i64)); }
                else { out.push_str(&format!("{:?}", x)); }
            }
            J::Str(t) => esc(t, out),
            J::Arr(a) => {
                if a.is_empty() { out.push_str("[]"); return; }
                // short numeric arrays stay on one line
                if a.iter().all(|v| matches!(v, J::Num(_) | J::Bool(_) | J::Null)) && a.len() <= 12 {
                    out.push('[');
                    for (i, v) in a.iter().enumerate() { if i > 0 { out.push_str(", "); } v.write(out, ind); }
                    out.push(']');
                    return;
                }
                out.push_str("[\n");
                for (i, v) in a.iter().enumerate() { out.push_str(&pad(ind + 1)); v.write(out, ind + 1); if i + 1 < a.len() { out.push(','); } out.push('\n'); }
                out.push_str(&pad(ind)); out.push(']');
            }
            J::Obj(o) => {
                if o.is_empty() { out.push_str("{}"); return; }
                out.push_str("{\n");
                for (i, (k, v)) in o.iter().enumerate() { out.push_str(&pad(ind + 1)); esc(k, out); out.push_str(": "); v.write(out, ind + 1); if i + 1 < o.len() { out.push(','); } out.push('\n'); }
                out.push_str(&pad(ind)); out.push('}');
            }
        }
    }

    pub fn parse(src: &str) -> Result<J, String> {
        let b: Vec<char> = src.chars().collect();
        let mut i = 0;
        let v = value(&b, &mut i)?;
        ws(&b, &mut i);
        if i != b.len() { return Err(format!("trailing characters at {}", i)); }
        Ok(v)
    }
}

impl From<f64> for J { fn from(x: f64) -> J { J::Num(x) } }
impl From<usize> for J { fn from(x: usize) -> J { J::Num(x as f64) } }
impl From<u64> for J { fn from(x: u64) -> J { J::Num(x as f64) } }
impl From<u32> for J { fn from(x: u32) -> J { J::Num(x as f64) } }
impl From<bool> for J { fn from(x: bool) -> J { J::Bool(x) } }
impl From<&str> for J { fn from(x: &str) -> J { J::Str(x.to_string()) } }
impl From<String> for J { fn from(x: String) -> J { J::Str(x) } }
impl<T: Into<J>> From<Vec<T>> for J { fn from(x: Vec<T>) -> J { J::Arr(x.into_iter().map(Into::into).collect()) } }
impl<T: Into<J>> From<Option<T>> for J { fn from(x: Option<T>) -> J { x.map_or(J::Null, Into::into) } }

fn esc(t: &str, out: &mut String) {
    out.push('"');
    for c in t.chars() {
        match c { '"' => out.push_str("\\\""), '\\' => out.push_str("\\\\"), '\n' => out.push_str("\\n"), '\t' => out.push_str("\\t"), '\r' => out.push_str("\\r"),
                  c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)), c => out.push(c) }
    }
    out.push('"');
}

fn ws(b: &[char], i: &mut usize) { while *i < b.len() && b[*i].is_whitespace() { *i += 1; } }

fn value(b: &[char], i: &mut usize) -> Result<J, String> {
    ws(b, i);
    let c = *b.get(*i).ok_or("unexpected end")?;
    match c {
        '{' => {
            *i += 1; let mut o = vec![];
            loop {
                ws(b, i);
                if b.get(*i) == Some(&'}') { *i += 1; break; }
                let k = match value(b, i)? { J::Str(s) => s, _ => return Err("object key must be a string".into()) };
                ws(b, i);
                if b.get(*i) != Some(&':') { return Err(format!("expected ':' at {}", i)); }
                *i += 1;
                o.push((k, value(b, i)?));
                ws(b, i);
                match b.get(*i) { Some(',') => *i += 1, Some('}') => { *i += 1; break; } _ => return Err(format!("expected ',' or '}}' at {}", i)) }
            }
            Ok(J::Obj(o))
        }
        '[' => {
            *i += 1; let mut a = vec![];
            loop {
                ws(b, i);
                if b.get(*i) == Some(&']') { *i += 1; break; }
                a.push(value(b, i)?);
                ws(b, i);
                match b.get(*i) { Some(',') => *i += 1, Some(']') => { *i += 1; break; } _ => return Err(format!("expected ',' or ']' at {}", i)) }
            }
            Ok(J::Arr(a))
        }
        '"' => {
            *i += 1; let mut s = String::new();
            loop {
                let c = *b.get(*i).ok_or("unterminated string")?; *i += 1;
                match c {
                    '"' => break,
                    '\\' => { let e = *b.get(*i).ok_or("bad escape")?; *i += 1; match e { 'n' => s.push('\n'), 't' => s.push('\t'), 'r' => s.push('\r'), '"' => s.push('"'), '\\' => s.push('\\'), '/' => s.push('/'),
                        'u' => { let h: String = b[*i..(*i + 4).min(b.len())].iter().collect(); *i += 4; s.push(char::from_u32(u32::from_str_radix(&h, 16).map_err(|e| e.to_string())?).unwrap_or('?')); }
                        _ => return Err("bad escape".into()) } }
                    c => s.push(c),
                }
            }
            Ok(J::Str(s))
        }
        't' if b[*i..].starts_with(&['t', 'r', 'u', 'e']) => { *i += 4; Ok(J::Bool(true)) }
        'f' if b[*i..].starts_with(&['f', 'a', 'l', 's', 'e']) => { *i += 5; Ok(J::Bool(false)) }
        'n' if b[*i..].starts_with(&['n', 'u', 'l', 'l']) => { *i += 4; Ok(J::Null) }
        _ => {
            let st = *i;
            while *i < b.len() && (b[*i].is_ascii_digit() || matches!(b[*i], '-' | '+' | '.' | 'e' | 'E')) { *i += 1; }
            let t: String = b[st..*i].iter().collect();
            t.parse::<f64>().map(J::Num).map_err(|_| format!("bad token '{}' at {}", t, st))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn roundtrip() {
        let j = J::obj().with("a", 1.5).with("b", "x\"y\n").with("c", vec![1.0, 2.0, 3.0]).with("d", J::obj().with("e", true).with("f", J::Null)).with("g", 12345678901u64);
        let t = j.pretty();
        let p = J::parse(&t).unwrap();
        assert_eq!(p.f("a"), Some(1.5));
        assert_eq!(p.s("b"), Some("x\"y\n"));
        assert_eq!(p.arr("c").unwrap().len(), 3);
        assert_eq!(p.get("d").unwrap().b("e"), Some(true));
        assert_eq!(p.f("g"), Some(12345678901.0));
        assert!(J::parse("{\"a\": }").is_err());
    }
    #[test]
    fn floats_keep_precision() {
        let x = 0.1 + 0.2;
        let t = J::obj().with("x", x).pretty();
        assert_eq!(J::parse(&t).unwrap().f("x"), Some(x));
    }
}
