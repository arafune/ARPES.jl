# Development Workflow

This project follows a specific branching strategy to comply with Julia General Registry requirements regarding submodules.

## Branching Strategy

- **`dev` (Default Branch):** All development happens here. This branch contains the `testdata` submodule for CI and local testing. **All Pull Requests should be targeted at `dev`.**
- **`main` (Deployment Branch):** This branch is an automated, "clean" mirror of `dev`. It **does not** contain submodules or `.gitmodules`. This branch is used for registration in the General Registry.

## Automated Sync (GitHub Actions)

We use a workflow (`.github/workflows/deploy-main.yml`) that triggers on every push to `dev`.

### What the workflow does

1. Checks out the `dev` branch.
2. Removes the `testdata` submodule and `.gitmodules` file.
3. Force-pushes the result to the `main` branch.

**Note:** Never commit directly to `main`. Any manual changes to `main` will be overwritten by the next sync from `dev`.
