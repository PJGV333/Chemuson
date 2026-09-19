# Design

M21 owns a small settings policy surface: a QSettings factory, a minimal
key/value protocol, boolean coercion, and typed naming/numbering preference
round-trips. It also owns the canonical `open_resource_path` implementation.

The GUI shell creates the store through `application_settings()` and passes it
to existing controllers. Main-window and view-controller wrappers retain their
historical method names and behavior while delegating normalization and
persistence to M21. The update controller receives the same structural settings
contract without importing QSettings directly.

`utils/resources.py` remains a compatibility shim. ChemName consumes the
canonical platform resource helper; old imports continue to resolve.
