# Burgers.jl

Julia package for solving Burgers' equation:
$$
\frac{\partial u}{\partial t}+ \frac{\partial}{\partial x}\left( \frac{u^2}{2} \right) = 0.
$$

# Installation

In Julia REPL, press `]` and type `add https://github.com/dmitry-kabanov/Burgers.jl`

# Development

For developing the package locally, press `]` in Julia REPL to switch to Package mode, and type:

    develop --local https://github.com/dmitry-kabanov/Burgers.jl
    
Type `Backspace` to return to Julia mode and type `using Burgers`.
Now every change to the package will be visible to the current REPL environment.
