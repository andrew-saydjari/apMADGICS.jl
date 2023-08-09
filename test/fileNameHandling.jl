@testset "fileNameHandling.jl" begin

    # Need to read through these to see if we think this is 
    # the optimal caching strategy

    rtest = (10, "lco25m", "000+04", "10155", "58656", "asVisit-dr17-10155-58656-151.fits", "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit/lco25m/000+04/10155/58656/asPlate-a-10155-58656.fits", 150)

    itest = (10, "apo25m0000010i", "180+60", "9667", "58126", "apVisit-dr17-9667-58126-002.fits", "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit/apo25m/180+60/9667/58126/apPlate-a-9667-58126.fits", 295)

    @test cache_skyname(rtest[2:end]) == "../local_cache/58656/sky_lco25m_000+04_10155_58656.jdat"
    @test cache_skyname(itest[2:end]) == "../local_cache/58126/sky_apo25m0000010i_180+60_9667_58126.jdat"

    @test cache_skynameSpec(rtest[2:end]) == "../local_cache/58656/450/skySpec_lco25m_000+04_10155_58656_150.jdat"
    @test cache_skynameSpec(itest[2:end]) == "../inject_local_cache/58126/295/skySpec_apo25m0000010i_180+60_9667_58126_295.jdat"

    @test cache_starname(rtest[2:end]) == "../local_cache/58656/450/star_lco25m_000+04_10155_58656_150.jdat"
    @test cache_starname(itest[2:end]) == "../inject_local_cache/58126/295/star_apo25m0000010i_180+60_9667_58126_295.jdat"

    # add frameFromPlate test here
    @test getFramesFromPlate("../data/asPlate-a-10155-58656.fits") == ["30940037", "30940038", "30940039", "30940040", "30940041", "30940042", "30940043", "30940044"]

    ## build_framepath tests here
    @test build_framepath(rtest[2],rtest[5],"30940037","a") == "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/exposures/apogee-s/58656/as1D-a-30940037.fits"

    #this does not exist, but just testing the n/s switch
    @test build_framepath(itest[2],itest[5],"30940037","a") =="/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/exposures/apogee-n/58126/ap1D-a-30940037.fits"

    vrname = build_visitpath(rtest[2:7])
    @test vrname == "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit/lco25m/000+04/10155/58656/asVisit-dr17-10155-58656-151.fits"
    viname = build_visitpath(itest[2:7])
    @test viname == 
    "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit/apo25m/180+60/9667/58126/apVisit-dr17-9667-58126-002.fits"

    ## visit2cframe test
    @test visit2cframe(vrname,rtest[2],"30940037","a") == "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit/lco25m/000+04/10155/58656/asCframe-a-30940037.fits"
    @test visit2cframe(viname,itest[2],"30940037","a") == "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit/apo25m/180+60/9667/58126/apCframe-a-30940037.fits" 

    @test plate2visit(rtest[7],rtest[8]) == "asVisit-dr17-10155-58656-151.fits"
    @test plate2visit(itest[7],itest[8]) == "apVisit-dr17-9667-58126-006.fits"
end