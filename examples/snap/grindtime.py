import subprocess

def run(ngrid1d):
    ngrid = ngrid1d**3
#    print(f"ngrid1d = {ngrid1d}")
#    print(f"nsteps = {nsteps}")
    print(f"ngrid = {ngrid}")

    exe = "../../build/lmp"
    arglist = ["-in", f"{infile}", "-var", "nsteps", f"{nsteps}", "-var", "ngrid", f"{ngrid1d}", "-var", "nrep", f"{nrep}", "-var", "a", f"{alat}", "-var", "rcutfac", f"{rcutfac}"]
    foo = subprocess.run([exe] + arglist, capture_output=True)

    #log = open("log.lammps",'r')
    #lines = log.readlines()

    for line in str(foo.stdout).split('\\n'):
        if "Loop" in line:
            words = line.split()
            looptime = float(words[3])
            ncomputes = int(words[8])

    grindtime = ncomputes*ngrid/looptime

    if ncomputes != nsteps:
        print("Error: ncomputes != nsteps")
        exit()

#    print(f"Loop time = {looptime} secs.")
#    print(f"Number of computes = {ncomputes}")
    print(f"Grind time = {grindtime} gridpoints/sec")

    return grindtime

nrep = 3
ngridmin = 2
ngridmax = 8
ngrid1d = ngridmin
nsteps = 10
infile = "in.grindtime"
alat = 4.0495
rcutfac = 4.67637

grindtime = run(ngrid1d)
for ngrid1d in range(ngridmin, ngridmax):
    grindtime = run(ngrid1d)
    
    
