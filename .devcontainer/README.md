# Dev Container for JuliARPES

This directory contains the configuration for a Development Container (Dev Container) that provides a fully configured Julia 1.12 development environment for the JuliARPES package.

## What is a Dev Container?

A Dev Container is a Docker container configured specifically for development. It includes all the necessary tools, extensions, and dependencies to work on this project without manual setup.

## Requirements

To use this Dev Container, you need:

1. **Docker Desktop** (or Docker Engine on Linux)
   - [Download Docker Desktop for Mac/Windows](https://www.docker.com/products/docker-desktop/)
   - On Linux, install Docker Engine via your package manager

2. **Visual Studio Code** with the Remote - Containers extension
   - [Download VS Code](https://code.visualstudio.com/)
   - Install the [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)

**OR**

3. **GitHub Codespaces** (no local installation required)
   - Available directly from the GitHub repository

## Quick Start

### Option 1: Using VS Code Locally

1. Make sure Docker is running on your machine
2. Open the repository in VS Code
3. When prompted, click **"Reopen in Container"**
   - Or use Command Palette (F1) → **"Dev Containers: Reopen in Container"**
4. Wait for the container to build and initialize (first time takes a few minutes)
5. Once ready, you'll have a fully configured Julia 1.12 environment!

### Option 2: Using GitHub Codespaces

1. Go to the repository on GitHub
2. Click the **Code** button → **Codespaces** tab
3. Click **"Create codespace on main"** (or your branch)
4. Wait for the environment to initialize
5. Start developing in your browser or connect from VS Code!

## What's Included

### Pre-installed Tools
- Julia 1.12 (compatible with the package requirements)
- Git for version control
- Zsh with Oh My Zsh for an enhanced terminal experience

### VS Code Extensions
- **Julia Language Support** - Syntax highlighting, IntelliSense, debugging
- **GitLens** - Enhanced Git integration
- **Code Spell Checker** - Catch typos in code and comments
- **Live Share** - Collaborative editing

### Automatic Setup
When the container is created, it automatically:
1. Instantiates the main package dependencies (`Pkg.instantiate()`)
2. Builds the package if necessary (`Pkg.build()`)
3. Sets up the test environment
4. Sets up the documentation environment

### Port Forwarding
The following ports are automatically forwarded:
- **8000-8001**: For HTTP servers
- **8050-8051**: For Pluto notebooks

## Working in the Container

### Running Tests
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

### Building Documentation
```bash
cd docs
julia --project make.jl
```

### Starting a Julia REPL
```bash
julia --project=.
```

### Using JuliaFormatter
The project uses the "blue" style for formatting:
```bash
julia --project=. -e 'using JuliaFormatter; format("src")'
```

## Customization

If you need to customize the Dev Container for your workflow:

1. Edit `.devcontainer/devcontainer.json` to:
   - Add more VS Code extensions
   - Modify VS Code settings
   - Add additional features or tools
   - Change post-create commands

2. Rebuild the container:
   - Command Palette (F1) → **"Dev Containers: Rebuild Container"**

## Troubleshooting

### Container fails to build
- Ensure Docker is running
- Check your internet connection (needed to download images)
- Try rebuilding: Command Palette → **"Dev Containers: Rebuild Container Without Cache"**

### Julia packages fail to install
- Check your internet connection
- The first installation might take longer due to precompilation
- Check the terminal output for specific error messages

### Performance issues
- Allocate more resources to Docker in Docker Desktop settings
- On Windows, ensure you're using WSL 2 backend for better performance

## Additional Resources

- [VS Code Dev Containers Documentation](https://code.visualstudio.com/docs/devcontainers/containers)
- [GitHub Codespaces Documentation](https://docs.github.com/en/codespaces)
- [Julia Documentation](https://docs.julialang.org/)
