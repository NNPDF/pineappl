use assert_cmd::Command;

const HELP_STR: &str = "Display a manpage for selected subcommands

Usage: pineappl help [SUBCOMMAND]...

Arguments:
  [SUBCOMMAND]...  Name of the (chain of) subcommand(s) to read the manpage of

Options:
  -h, --help  Print help
";

const HELP_CONVOLUTE_STR: &str = r#".ie \n(.g .ds Aq \(aq
.el .ds Aq '
.TH convolute 1  "convolute " 
.SH NAME
convolute \- Convolutes a PineAPPL grid with a PDF set
.SH SYNOPSIS
\fBconvolute\fR [\fB\-b\fR|\fB\-\-bins\fR] [\fB\-i\fR|\fB\-\-integrated\fR] [\fB\-o\fR|\fB\-\-orders\fR] [\fB\-\-digits\-abs\fR] [\fB\-\-digits\-rel\fR] [\fB\-h\fR|\fB\-\-help\fR] <\fIINPUT\fR> <\fIPDFSETS\fR> 
.SH DESCRIPTION
Convolutes a PineAPPL grid with a PDF set
.SH OPTIONS
.TP
\fB\-b\fR, \fB\-\-bins\fR=\fIBINS\fR
Selects a subset of bins
.TP
\fB\-i\fR, \fB\-\-integrated\fR=\fIINTEGRATED\fR
Show integrated numbers (without bin widths) instead of differential ones
.TP
\fB\-o\fR, \fB\-\-orders\fR=\fIORDERS\fR
Select orders manually
.TP
\fB\-\-digits\-abs\fR=\fIABS\fR [default: 7]
Set the number of fractional digits shown for absolute numbers
.TP
\fB\-\-digits\-rel\fR=\fIREL\fR [default: 2]
Set the number of fractional digits shown for relative numbers
.TP
\fB\-h\fR, \fB\-\-help\fR
Print help
.TP
<\fIINPUT\fR>
Path of the input grid
.TP
<\fIPDFSETS\fR>
LHAPDF id(s) or name of the PDF set(s)
"#;

#[test]
fn help() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["help", "--help"])
        .assert()
        .success()
        .stdout(HELP_STR);
}

#[test]
fn help_convolute() {
    Command::cargo_bin("pineappl")
        .unwrap()
        .args(["help", "convolute"])
        .assert()
        .success()
        .stdout(HELP_CONVOLUTE_STR);
}
