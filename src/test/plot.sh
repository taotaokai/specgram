#!/bin/bash

name=seis

../../bin/specgram in=${name}.sac out=${name}.nc fmin=0.1 fmax=50 nf=100 dt=0.1

gmt set \
  TIME_INTERVAL_FRACTION 0.05 \
  FORMAT_DATE_MAP yyyy-mm-dd \
  FORMAT_CLOCK_MAP hh:mm \
  FONT_ANNOT_PRIMARY +10p \
  FONT_ANNOT_SECONDARY +12p

cat <<EOF > log10_annots.txt
-2         a   10@+-2
-1.69897   f
-1.52288   f
-1.39794   f
-1.30103   f
-1.22185   f
-1.1549    f
-1.09691   f
-1.04576   f
-1         a   10@+-1
-0.69897   f
-0.522879  f
-0.39794   f
-0.30103   f
-0.221849  f
-0.154902  f
-0.09691   f
-0.0457575 f
0          a   10@+0
0.30103    f
0.477121   f
0.60206    f
0.69897    f
0.778151   f
0.845098   f
0.90309    f
0.954243   f
1          a   10@+1
1.30103    f
1.47712    f
1.60206    f
1.69897    f
1.77815    f
1.8451     f
1.90309    f
1.95424    f
2          a   10@+2
2.30103    f
2.47712    f
2.60206    f
2.69897    f
2.77815    f
2.8451     f
2.90309    f
2.95424    f
3          a   10@+3
EOF



gmt grdreformat ${name}.nc?stf_abs ${name}.grd

gmt grd2cpt ${name}.grd -Cjet -Z > ${name}.cpt

gmt grdimage ${name}.grd \
  -JX10i/6i -C${name}.cpt \
  -BWSen+t"time-frequency power plot" \
  -By+l"Frequency (Hz)" -Byc"log10_annots.txt" \
  -Bxa5Sf1S -K > ${name}.ps

gmt psscale -D5i/-0.5i/3i/0.1ih -Bxa1f0.5 -By+l"magnitude" -C${name}.cpt -O >> ${name}.ps

#gmt grdimage ${name}.grd \
#  -JX10i/6i -Cjet \
#  -R2014-05-18T00/2014-05-18T03/-1/2.48 \
#  -BWSen+t"time-frequency power plot" \
#  -By+l"Frequency (Hz)" -Byc"log10_annots.txt" \
#  -Bpxa30Mf10m -Bsxa1D > ${name}.ps

open ${name}.ps
