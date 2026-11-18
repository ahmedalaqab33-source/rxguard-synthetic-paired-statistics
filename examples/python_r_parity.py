"""Compare two implementations on eight invented binary rows and a shared schedule."""
from pathlib import Path
import argparse,csv,json,subprocess,sys,math
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from paired_stats_v1 import invented_rows,resampling_plan,cluster_bootstrap,paired_summary
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--output',type=Path,required=True)
 p.add_argument('--rscript',required=True)
 a=p.parse_args();out=a.output.resolve()
 if out==ROOT or ROOT in out.parents:raise SystemExit('Choose a new output outside the repository')
 out.mkdir(parents=True,exist_ok=False)
 rows=invented_rows();plan=resampling_plan(rows,draws=2000,synthetic=True)
 with (out/'synthetic_resampling_plan.csv').open('x',newline='',encoding='utf-8') as f:
  w=csv.writer(f);w.writerow(['scope','replicate','draw','cluster_id'])
  for i,replicate in enumerate(plan,1):
   for j,cluster in enumerate(replicate,1):w.writerow(['SYNTHETIC',i,j,cluster])
 python=cluster_bootstrap(rows,plan,synthetic=True)
 proc=subprocess.run([a.rscript,'--vanilla',str(ROOT/'paired_stats_v1.R'),str(out)],
     capture_output=True,text=True,timeout=120,check=False)
 (out/'r_stdout.txt').write_text(proc.stdout,encoding='utf-8')
 (out/'r_stderr.txt').write_text(proc.stderr,encoding='utf-8')
 if proc.returncode:raise SystemExit('R comparison failed; inspect external logs')
 with (out/'r_bootstrap_values.csv').open(newline='') as f:
  values=[float(r['value']) for r in csv.DictReader(f)]
 r=json.loads((out/'r_mechanics.json').read_text())
 assert len(values)==len(python['values'])==2000
 assert all(math.isclose(x,y,rel_tol=0,abs_tol=1e-12) for x,y in zip(values,python['values']))
 assert math.isclose(python['lower'],r['bootstrap_lower'],abs_tol=1e-12)
 assert math.isclose(python['upper'],r['bootstrap_upper'],abs_tol=1e-12)
 result={'scope':'INVENTED_ARITHMETIC_ONLY','status':'PASS','replicates_compared':2000,
         'clinical_validation':False,'positive_stratum':paired_summary(rows,synthetic=True)}
 (out/'PARITY_SUMMARY.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print(json.dumps(result))
if __name__=='__main__':main()
