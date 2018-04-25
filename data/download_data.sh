for league in E D I SP F
do
    for i in $(seq 93 111)
    do
        if [ "$league" == "E" ]; then 
            div=0
        else
	    div=1
	fi
        m=$(printf %02d $((i % 100)))
        n=$(printf %02d $(((i+1)%100)))
        echo http://www.football-data.co.uk/mmz4281/"$m""$n"/"$league""$div".csv
        wget http://www.football-data.co.uk/mmz4281/"$m""$n"/"$league""$div".csv -O "$league"_"$m""$n".csv
        sleep 1s
    done
done
