# MacroDR CMD — Command Surfaces

- Purpose: stable, script-callable command APIs and small registry builders.
- What goes here: headers only; argument schemas and pre/postcondition stubs. No heavy logic.
- Implementations live in `src/core/` and may wrap legacy while migrating.
- Registry aggregation happens in `include/macrodr/cli/command_manager.h`.
- See `docs/architecture/modules.md` for boundaries and allowed dependencies.

