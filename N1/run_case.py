#!/usr/bin/env python3
"""usage: run_case.py <case dir>  -- run the R6 solver there (single thread) and record wall time, CPU time and peak memory in timing.json"""
import subprocess, sys, time, resource, json
BIN='/home/claude/work/r4/corrected/target/release/DES_thermal_simulation'
def run_case(d):
    t=time.time()
    r0=resource.getrusage(resource.RUSAGE_CHILDREN)
    with open(d+'/run.log','w') as lg: r=subprocess.run([BIN],cwd=d,stdout=lg,stderr=subprocess.STDOUT)
    ru=resource.getrusage(resource.RUSAGE_CHILDREN)
    info=dict(returncode=r.returncode,wall_s=time.time()-t,cpu_s=(ru.ru_utime+ru.ru_stime)-(r0.ru_utime+r0.ru_stime),peak_rss_MB=ru.ru_maxrss/1024)
    json.dump(info,open(d+'/timing.json','w')); return info
if __name__=='__main__': sys.exit(run_case(sys.argv[1])['returncode'])
