using TOML

if Sys.WORD_SIZE == 64 && (Sys.islinux() || Sys.isapple()) && isempty(VERSION.prerelease)
    root = dirname(@__DIR__)
    path = joinpath(root, "test", "Project.toml")
    project = TOML.parsefile(path)
    project["deps"]["Reactant"] = TOML.parsefile(joinpath(root, "Project.toml"))["weakdeps"]["Reactant"]
    open(io -> TOML.print(io, project; sorted = true), path, "w")
end
