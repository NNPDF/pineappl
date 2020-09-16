use prettytable::format::{FormatBuilder, LinePosition, LineSeparator};
use prettytable::Table;

pub(crate) fn create_table() -> Table {
    let mut table = Table::new();
    table.set_format(
        FormatBuilder::new()
            .column_separator(' ')
            .separator(LinePosition::Title, LineSeparator::new('-', '+', ' ', ' '))
            .build(),
    );
    table
}
