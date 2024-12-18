# Настройки терминала
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'

# Настройка вывода первого графика
set output 'velocity_vs_time.png'
set title 'Velocity vs Time'
set xlabel 'Time (s)'
set ylabel 'Velocity (m/s)'
set grid
plot 'data.txt' using 2:3 with lines title 'Velocity'

# Настройка вывода второго графика
set output 'distance_vs_time.png'
set title 'Distance vs Time'
set xlabel 'Time (s)'
set ylabel 'Distance (m)'
set grid
plot 'data.txt' using 2:4 with lines title 'Distance'
