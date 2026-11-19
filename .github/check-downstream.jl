using Pkg

# Resolve the candidate and its source-pinned dependencies together in
# the downstream environment before building or running its tests.
dev = [Pkg.PackageSpec(path = dirname(Base.active_project()))]
for dep in values(Pkg.dependencies())
    if dep.is_direct_dep && (dep.is_tracking_path || dep.is_tracking_repo)
        push!(dev, Pkg.PackageSpec(path = dep.source))
    end
end
Pkg.activate(pwd())
Pkg.develop(dev)
Pkg.instantiate()
Pkg.build()
Pkg.test()
