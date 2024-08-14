Implementation
==============

The wrapper is built using PyO3 to interface Rust to Python.

- The original Rust library that actually provides all the functionality. Any true code should be placed here as this is shared by all interfaces.
- The PyO3 wrapper objects written in Rust that define which methods are exposed to Python and that takes care of the type conversion from Python-understandable types to Rust types.
- The PyO3 wrapper objects in Python which are an exact mirror of their Rust equivalent. This translation is provided by PyO3 and note that it also preserves the associated documentation.
