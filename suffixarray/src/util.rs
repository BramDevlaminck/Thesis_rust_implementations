use std::time::{SystemTime, SystemTimeError, UNIX_EPOCH};

use crate::searcher::Searcher;

/// Returns the current UNIX time in milliseconds
pub fn get_time_ms() -> Result<f64, SystemTimeError> {
    Ok(SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos() as f64 * 1e-6)
}

/// Returns how long executing `f` took in ms
pub fn time_execution(
    searcher: &mut Searcher,
    f: &dyn Fn(&mut Searcher) -> bool,
) -> Result<(bool, f64), SystemTimeError> {
    let start_ms = get_time_ms()?;
    let found = f(searcher);
    let end_ms = get_time_ms()?;
    Ok((found, end_ms - start_ms))
}
