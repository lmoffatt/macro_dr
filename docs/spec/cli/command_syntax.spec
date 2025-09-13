Here is a draft for `spec/cli/command_syntax.spec`, capturing the current and planned command-line structure of MacroIR:

---

````spec
# MacroIR Command Syntax Specification

Version: 0.1  
Status: DRAFT  
Last updated: 2025-07-17  
Scope: CLI interface only — excludes internal DSL syntax

---

## Top-level Structure

macroir <command> [<target>] [--<flag>=<value>]...

Where:
- `<command>` is a required verb describing the high-level operation.
- `<target>` is one or more arguments specific to the command.
- `--<flag>=<value>` are optional global or command-specific modifiers.

---

## Supported Commands

### 1. run

```sh
macroir run <file1.macroir> [file2.macroir ...] [--flags]
````

* Executes one or more `.macroir` DSL scripts sequentially.
* Flags may configure execution (e.g., verbosity, export mode).
* If `<command>` is omitted, `run` is assumed.

### 2. eval

```sh
macroir eval "<expression>"
```

* Evaluates a single DSL expression passed as a string.
* Returns output to stdout (unless `--quiet` is set).

### 3. compile

```sh
macroir compile <file.macroir>
```

* Parses and type-checks a DSL script without executing it.
* Useful for validation, static analysis, or ahead-of-time compilation.

### 4. describe

```sh
macroir describe             # Lists all known DSL functions
macroir describe <function> # Details about a single function
```

* Used for introspection and developer documentation.

### 5. help

```sh
macroir help
```

* Displays this help text or command-specific guidance.

---

## Flags (global or per-command)

| Flag            | Type    | Applies to | Description                             |
| --------------- | ------- | ---------- | --------------------------------------- |
| `--json`        | boolean | all        | Emit results in JSON format             |
| `--dry-run`     | boolean | run, eval  | Parse and compile only, do not execute  |
| `--log-level`   | string  | all        | One of `error`, `warn`, `info`, `debug` |
| `--output=FILE` | string  | run, eval  | Redirect main result to a file          |

---

## Aliases and Default Behavior

* If no `<command>` is provided and the first argument is a `.macroir` file, it is treated as:

  ```sh
  macroir run script.macroir
  ```

---

## Reserved Extensions (for command target inference)

| Extension       | Meaning                       |
| --------------- | ----------------------------- |
| `.macroir`      | DSL script                    |
| `.macroir.json` | Environment/config definition |
| `.macroir.csv`  | Structured output file        |
| `.macroir.log`  | Execution trace               |

---

## Examples

```sh
macroir run model.macroir
macroir eval "simulate(model, parameters)"
macroir compile scripts/test.macroir
macroir describe simulate
macroir help
```

---

## Future Extensions (planned)

* `check`: pre-run validation of input/output structure
* `format`: canonicalize DSL script formatting
* `init`: generate a boilerplate environment or script

```

---

Let me know if you'd like to auto-generate CLI help output from this spec!
```

