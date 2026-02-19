using Test
using ARPES

@testset "_parse_waves_declaration basic test" begin
    decl = "WAVES/D/N=5 testwave"
    result = ARPES.IO.Format._parse_waves_declaration(decl)
    @test result[:name] == "testwave"
    @test result[:n] == 5
    @test result[:type] == :double_precision
end

@testset "_parse_waves_declaration single precision, quoted name, shape" begin
    decl = "WAVES/S/N=(200,2401) 'GrIr111_1'"
    result = ARPES.IO.Format._parse_waves_declaration(decl)
    @test result[:name] == "GrIr111_1"
    @test result[:n] == 200
    @test result[:type] == :single_precision
    @test result[:shape] == (200, 2401)
end

@testset "_parse_waves_declaration default type, unquoted name" begin
    decl = "WAVES/N=10 mywave"
    result = ARPES.IO.Format._parse_waves_declaration(decl)
    @test result[:name] == "mywave"
    @test result[:n] == 10
    @test result[:type] == :double_precision
end


#setscale line format:
@testset "_parse_setscale! basic test" begin
    line = "X SetScale/I x, -11.44, 12.44, \"deg (theta_y)\", 'GrIr111_1'"
    waves_info = ARPES.IO.Format._parse_setscale!(line)
    @test haskey(waves_info, :x_scale)
    xscale = waves_info[:x_scale]
    @test xscale[:min] ≈ -11.44
    @test xscale[:max] ≈ 12.44
    @test xscale[:unit] == "deg"
    @test xscale[:label] == "theta_y"
end

@testset "_parse_setscale! no label" begin
    line = "X SetScale/I y, 0, 100, \"eV\", 'wave1'"
    waves_info = ARPES.IO.Format._parse_setscale!(line)
    @test haskey(waves_info, :y_scale)
    yscale = waves_info[:y_scale]
    @test yscale[:min] == 0
    @test yscale[:max] == 100
    @test yscale[:unit] == "eV"
    @test isnothing(yscale[:label])
end

@testset "_parse_setscale! " begin
    line = "X SetScale/I x, -10.6875, 10.6875, \"deg (Non-Energy Channel)\", 'Region1_1'"
    waves_info = ARPES.IO.Format._parse_setscale!(line)
    @test haskey(waves_info, :x_scale)
    xscale = waves_info[:x_scale]
    @test xscale[:min] ≈ -10.6875
    @test xscale[:max] ≈ 10.6875
    @test xscale[:unit] == "deg"
    @test xscale[:label] == "non_energy_channel"
end

@testset "_parse_setscale! basic test perpoints" begin
    line = "X SetScale/P y, 5.2, 0.002, \"eV\", 'GrIr111_1'"
    waves_info = ARPES.IO.Format._parse_setscale!(line)
    @test haskey(waves_info, :y_scale)
    yscale = waves_info[:y_scale]
    @test yscale[:start] ≈ 5.2
    @test yscale[:step] ≈ 0.002
    @test yscale[:unit] == "eV"
    @test isnothing(yscale[:label])
end

# comment line:
@testset"_parse_comment_line created_date" begin
    line = "Created Date (UTC): 2025-Dec-25 10:04:00.189347"
    result = ARPES.IO.Format._parse_comment_line(line)
    @test result[:created_date] == "2025-Dec-25 10:04:00.189347"
end

@testset"_parse_comment_line created_by" begin
    line = "Created by: SpecsLab Prodigy, Version 4.123.1-r123826 "
    result = ARPES.IO.Format._parse_comment_line(line)
    @test result[:created_by] == "SpecsLab Prodigy, Version 4.123.1-r123826"
end

@testset"_parse_comment_line scan_mode" begin
    line = "Scan Mode         = Fixed Analyzer Transmission"
    result = ARPES.IO.Format._parse_comment_line(line)
    @test result[:scan_mode] == "Fixed Analyzer Transmission"
end

@testset"_parse_comment_line number_of_scan" begin
    line = "Number of Scans   = 16"
    result = ARPES.IO.Format._parse_comment_line(line)
    @test result[:number_of_scans] == 16
end

@testset"_parse_comment_line photon energy" begin
    line = "Excitation Energy = 4.708"
    result = ARPES.IO.Format._parse_comment_line(line)
    @test result[:excitation_energy] ≈ 4.708
end

@testset"_parse_comment_line not include" begin
    line = "Acquisition Parameters:"
    result = ARPES.IO.Format._parse_comment_line(line)
    @test isnothing(result)
end

@testset"_parse_comment_line comment" begin
    line = "User Comment      = beta:0;Temperature:RT;X:13.5;Y:21.55;Z:+00346000;theta:-00270000;position:187.4655;UV(P);IR(P);+w;P:113mW;+3w;P:6mW;"
    result = ARPES.IO.Format._parse_comment_line(line)
    @test result[:user_comment] ==
          "beta:0;Temperature:RT;X:13.5;Y:21.55;Z:+00346000;theta:-00270000;position:187.4655;UV(P);IR(P);+w;P:113mW;+3w;P:6mW;"
    @test result[:beta] == 0
    @test result[:temperature] == "RT"
    @test result[:x] ≈ 13.5
    @test result[:z] ≈ 346000
    @test result[:position] ≈ 187.4655
end

@testset"_parse_comment_line_comment" begin
    line = "beta:0;Temperature:RT;X:13.5;Y:21.55;Z:+00346000;theta:-00270000;position:187.4655;UV(P);IR(P);+w;P:113mW;+3w;P:6mW;"
    result = ARPES.IO.Format._parse_comment_line_comment(line)
    @test result[:beta] == 0
    @test result[:temperature] == "RT"
    @test result[:x] ≈ 13.5
    @test result[:z] ≈ +00346000
    @test result[:position] ≈ 187.4655
end

@testset"_parse_wave_data" begin
    line1 = "67.8969 67.8969 71.7938 75.7938 77.8969 83.6906 84 85.8969 86 86 87.8969"
    line2 = ""
    line3 = "54.8032 54.8032 57.6065 60 60 63.213 69.6065 71.1968"
    data_lines = String[]
    push!(data_lines, line1)
    push!(data_lines, line2)
    push!(data_lines, line3)
    result = ARPES.IO.Format._parse_wave_data(data_lines)
    @test result[1] ≈ 67.8969
    @test length(result) == 19
    @test result[19] ≈ 71.1968
end
