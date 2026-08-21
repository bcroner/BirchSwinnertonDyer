from complete2 import complete_descent
print("  COMPLETE 2-DESCENT SWEEP  y^2 = x^3 - n^2 x")
print("   n | #Sel2 | dimSel2 | rank>= | rank<= | RANK | dim Sha[2] | unk")
print("  " + "-"*70)
res={}
for n in range(1,61):
    S,c,sel,img,unk = complete_descent(0,n,-n)
    ds=len(sel); dim=ds.bit_length()-1
    di=len(img); dimi=di.bit_length()-1
    lo=max(dimi-2,0); hi=dim-2
    pw = ds and not (ds&(ds-1))
    exact = (lo==hi)
    sha = (hi-lo) if exact else None
    res[n]=(ds,dim,lo,hi,len(unk),pw)
    tag = str(lo) if exact else f"[{lo},{hi}]"
    shas = "0" if exact else f"<={hi-lo}"
    print(f"  {n:>3} | {ds:>5} | {dim:>7} | {lo:>6} | {hi:>6} | {tag:>6} | {shas:>10} | {len(unk):>3}"
          f"{'' if pw else '   NOT-2^k'}")
ex=[n for n in res if res[n][2]==res[n][3]]
gp=[n for n in res if res[n][2]!=res[n][3]]
bad=[n for n in res if not res[n][5]]
un=[n for n in res if res[n][4]]
print(f"\n  exact ranks : {len(ex)}/60")
print(f"  still a gap : {gp}")
print(f"  non-power-of-2 (bug indicator): {bad}")
print(f"  inconclusive local tests      : {un}")
