# run for a set of instances
ls -1 instances/*_{00008,00016,00032}_*.json | xargs -n1 -I{} julia --project src/main.jl {} 0.6 \;
# match with the IP and LP objective values
tail -1 -q instances/*.txt | sort | join - best_ip_lp_obj.txt >results.txt
