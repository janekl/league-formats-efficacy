for league in D E F I SC SP
do
    for i in $(seq 2 16)
    do
        if [ "$league" == "E" ] || [ "$league" == "SC" ]; then 
            div=0
        else
            div=1
        fi
        m=$(printf %02d $i)
        n=$(printf %02d $((i+1)))
        echo http://www.football-data.co.uk/mmz4281/"$m""$n"/"$league""$div".csv
        wget http://www.football-data.co.uk/mmz4281/"$m""$n"/"$league""$div".csv -O "$league"_"$m""$n".csv
        sleep 1s
    done
done
wget http://www.football-data.co.uk/new/POL.csv
