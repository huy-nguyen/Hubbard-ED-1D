function test_suite = test_equal_time_gf
initTestSuite;
end

function inp = setup

inp.in_t = 1;
inp.in_U = 4;
inp.method = 'short_tau';
inp.commit_number = 'testtesttest';
inp.sector = 'up';
inp.NUM_CORES = 4;
end

function test_5_sites(inp)
in_noOfSites = 5;
in_noOfUp = 3;
in_noOfDn = 3;
output_file =  equalTimeGF( inp.in_t, inp.in_U, in_noOfSites, in_noOfUp, in_noOfDn, inp.sector, inp.method, inp.commit_number, inp.NUM_CORES );
load(output_file, '-mat', 'spinUpGreenFunction');
assertElementsAlmostEqual(spinUpGreenFunction, ...
    [0.4 -0.28245 0.10829 0.10829 -0.28245],...
       'relative', 0.0001);
end

function test_8_sites(inp)
in_noOfSites = 8;
in_noOfUp = 5;
in_noOfDn = 5;
output_file =  equalTimeGF( inp.in_t, inp.in_U, in_noOfSites, in_noOfUp, in_noOfDn, inp.sector, inp.method, inp.commit_number, inp.NUM_CORES );
load(output_file, '-mat', 'spinUpGreenFunction');
assertElementsAlmostEqual(spinUpGreenFunction, ...
    [0.375 -0.264695877699632 0.107783367516430 0.04369396434612655 -0.101808620724713 0.04369396432796521 0.107783367525266 -0.264695877690294],...
       'relative', 0.0001);
end