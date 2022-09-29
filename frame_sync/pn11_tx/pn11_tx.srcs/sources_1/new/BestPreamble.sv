`timescale 1ns / 1ps
`default_nettype none
//////////////////////////////////////////////////////////////////////////////////
// Company: ICE lab
// Engineer: Preston Walker
// 
// Create Date: 09/29/2022 11:18 AM
// Design Name: 
// Module Name: Best Preamble
// Project Name: Generate PN11
// Target Devices: 
// Description: 
// 
// Dependencies: 
// 
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module BestPreamble(
    input wire logic clk,
    input wire logic start,
    output logic finished,
    output logic out_bit
);

logic[6:0] counter128;
logic[127:0] preamble_bits = 128'b10011111011110001100110110111010011011110110100010100100101101011010101011010011110010110111011001010010111000101000011011100101;


always_ff @(posedge clk)
begin 
    // if the start goes high, reset the counters
    if (start) 
    begin 
        counter128 <= 127;
        finished <= 0;
        out_bit <= preamble_bits[counter128];
    end

    // output the bits in the preamble one at a time
    else
    begin
        out_bit <= preamble_bits[counter128-1];
        counter128 <= counter128 - 1;
            
        if (counter128 == 1)
        begin
            finished <= 1;
        end
    end
end


endmodule