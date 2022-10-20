add_force CLK100MHZ {0 0} {1 5ns} -repeat_every 10ns

add_force btnC 0 
add_force btnU 1

run 95 ns

add_force btnC 1
add_force btnU 0

run 50000 ns
