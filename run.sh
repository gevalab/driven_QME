name='no_field'

for (( val=0; val<=10; ++val )); do

# calculate parameters
w01=0.0294
k=0.2
w=$( echo "$w01*$val*$k" | bc)
epsilon=0.0
echo "number $val created:"
echo "w = $w"	
echo "epsilon = $epsilon"

# generate folders
mkdir $name$val
cd $name$val

# create input file
cat > input.txt <<!
projector           full
time_convolution    tc
spectral_density    Ohmic_2
steady_state        no
#############
w       $w
w01     0.0294
epsilon $epsilon
Lambda  1.102e-3
N       101
Omega   3.507e-4
y0      3.116e2
wc      3.507e-4
eta     1.066e6
beta    1.053e3
M       4.529e8
N_point 400
#############
T       2.067e4
dt      10.0
#############
sigma
1   0
0   0

end

!

# run the code
/home/zongweih/2022/driven_QME-main/driven_QME &
cd ../
done
