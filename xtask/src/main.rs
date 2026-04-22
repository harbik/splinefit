use std::process::{Command, exit};

fn main() {
    match std::env::args().nth(1).as_deref() {
        Some("test") => task_test(),
        Some("doc")  => task_doc(),
        _ => {
            eprintln!("Usage: cargo xtask <task>\n");
            eprintln!("Tasks:");
            eprintln!("  test    Run all tests");
            eprintln!("  doc     Regenerate README.md from lib.rs docs, then run doc-tests");
            exit(1);
        }
    }
}

/// Run all tests (unit, integration, and doc-tests).
fn task_test() {
    run("cargo", &["test"]);
}

/// Regenerate README.md via `cargo rdme`, then run doc-tests to verify the examples.
fn task_doc() {
    run("cargo", &["rdme", "--force"]);
    run("cargo", &["test", "--doc"]);
}

fn run(cmd: &str, args: &[&str]) {
    let status = Command::new(cmd)
        .args(args)
        .status()
        .unwrap_or_else(|e| {
            eprintln!("error: could not run `{} {}`: {e}", cmd, args.join(" "));
            exit(1);
        });
    if !status.success() {
        exit(status.code().unwrap_or(1));
    }
}
