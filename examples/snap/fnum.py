import lammps

infile_sna = "in.fnum_sna"
infile_gaussian = "in.fnum_gaussian"

fdeltalist = [1.0e-2,1.0e-3,1.0e-4,1.0e-5,1.0e-6,1.0e-7,1.0e-8,1.0e-9,1.0e-10]

print("Fdelta RMSE(SNAP) RMSE(Gaussian)")
for fdelta in fdeltalist:
    cmdlist = ["-screen","none","-var","fdelta",f'{fdelta}']
    lmp = lammps.lammps(cmdargs = cmdlist)
    lmp.file(infile_sna)
    faverrsq = lmp.extract_compute("faverrsq", lammps.LMP_STYLE_GLOBAL, lammps.LMP_TYPE_SCALAR)
    rmse_sna = faverrsq**0.5
    lmp.command("clear")
    lmp.file(infile_gaussian)
    faverrsq = lmp.extract_compute("faverrsq", lammps.LMP_STYLE_GLOBAL, lammps.LMP_TYPE_SCALAR)
    rmse_gaussian = faverrsq**0.5
    print(f"{fdelta} {rmse_sna} {rmse_gaussian}")
