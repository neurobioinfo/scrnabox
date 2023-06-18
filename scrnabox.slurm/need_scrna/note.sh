bash aa.sh -e aa -o bb.log
bash aa.sh -e aa -o bb.log 2>&1|tee -a as.log
bash aa.sh --help

bash a1.sh --help

bash a1.sh -h 
bash /Users/sam/aa.sh -h 


bash aa0.sh -e aa -o bb.log 2>&1|tee -a as2.log

&>
bash aa0.sh -e aa -o bb.log &> as3.log

cat as3.log

echo "hi" | grep "h"