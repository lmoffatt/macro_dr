


# 📁 `spec/` — MacroIR Specification Directory

This directory contains all formal specifications, documentation, and schemas that define the behavior and structure of the MacroIR system. It is organized into subdirectories by concern: CLI interface, DSL syntax, type rules, runtime contracts, and data schemas.

---

## 📦 Structure and Naming Conventions

| File/Folder                    | Purpose                                         | Format           |
|-------------------------------|--------------------------------------------------|------------------|
| `cli/command_syntax.spec`     | CLI verb interface and argument syntax          | `.spec` (plain)  |
| `dsl/grammar.bnf`             | Core grammar for the MacroIR DSL                | `.bnf`           |
| `dsl/type_rules.spec`         | Type deduction rules for DSL functions          | `.spec`          |
| `dsl/functions.md`            | User-facing documentation of DSL functions      | `.md`            |
| `schema/environment.schema.json` | JSON Schema for `.macroir.json` config files | `.schema.json`   |
| `postconditions/postconditions.dsl` | DSL for runtime validation conditions     | `.dsl` (custom)  |

---

## 🧠 Naming Rules

- `.spec` files define **human-readable** specifications of syntax, typing, or runtime contracts.
- `.bnf` and `.ebnf` files define **machine-readable** grammars for parsing the DSL.
- `.schema.json` defines **structured validation** for JSON-based configuration files.
- `.md` is used for **developer-friendly documentation**, like internal function catalogs.
- `.dsl` denotes **embedded or specialized DSL fragments** (e.g., for postcondition logic).

---

## 📌 Summary

The `spec/` directory ensures the MacroIR system remains:
- Consistent across CLI, DSL, and runtime behavior
- Extendable with formal schemas and compositional logic
- Documented for both users and contributors

```

---

Would you like me to generate a minimal `functions.md` or populate `command_syntax.spec` with the earlier CLI description?
