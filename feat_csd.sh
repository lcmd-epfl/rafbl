OUT="ligs/csd_pool.csv"

echo "name,c_smiles,nbo1,nbo2,onbo1,onbo2,cnbo1,cnbo2,angle,homo,lumo,dipole,na,na_an,nb,bpa,estrada,wiener,global_eff,balaban,hosoya,zagreb1,zagreb2,global_simple,crest_flex,kier_a,kier_b,kier_al,0k,1k,2k,3k,1ka,2ka,3ka,1kb,2kb,3kb,1kal,2kal,3kal,k_phi,k_phia,k_xia,k_phib,k_xib,k_phial,k_xial,redu,0chi,1chi,2chi,3chi,4chi,5chi,0chiv,1chiv,2chiv,3chiv,4chiv,5chiv,T,Tnbo1,Tnbo2,Tenbo1,Tenbo2,t_vol,t_sur,t_ova,t_vol_nh,t_sur_nh,t_ova_nh,b_vol,b_vol_nh,bonh1,bonh2,bonh3,bonh4,bonh5,bonh6,bonh7,bonh8,o1,o2,o3,o4,o5,o6,o7,o8,onh1,onh2,onh3,onh4,onh5,onh6,onh7,onh8" > $OUT

for f in ligs/csd_pool/*.log
do python main_feats.py $f 1 >> $OUT
done
