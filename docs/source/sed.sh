#Fix Image links
sed -i 's/\!\[alt tag\](https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/_img\/day1_morning\//![alt tag](/g' day1_morning.md
sed -i 's/\!\[alt tag\](https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/_img\/day2_morning\//![alt tag](/g' day2_morning.md
sed -i 's/\!\[alt tag\](https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/_img\/day3_morning\//![alt tag](/g' day3_morning.md
sed -i 's/\!\[alt tag\](https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/_img\/day1_afternoon\//![alt tag](/g' day1_afternoon.md
sed -i 's/\!\[alt tag\](https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/_img\/day2_afternoon\//![alt tag](/g' day2_afternoon.md
sed -i 's/\!\[alt tag\](https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/_img\/day3_afternoon\//![alt tag](/g' day3_afternoon.md
sed -i 's/\!\[alt tag\](https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/_img\/day1_after\//![alt tag](/g' day1_afternoon.md

#Fix Home links
sed -i 's/https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/README.md/index.html/g' *.md

#Fix Back to top link
sed -i 's/https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/day1_morning\/README.md/day1_morning.html/g' *.md
sed -i 's/https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/day2_morning\/README.md/day2_morning.html/g' *.md
sed -i 's/https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/day3_morning\/README.md/day3_morning.html/g' *.md
sed -i 's/https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/day1_afternoon\/README.md/day1_afternoon.html/g' *.md
sed -i 's/https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/day2_afternoon\/README.md/day2_afternoon.html/g' *.md
sed -i 's/https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/day3_afternoon\/README.md/day3_afternoon.html/g' *.md

#Fix online resources
sed -i 's/https:\/\/github.com\/alipirani88\/Comparative_Genomics\/blob\/master\/online_resources\/README.md/online_resources.html/g' *.md




