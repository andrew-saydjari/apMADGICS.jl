@testset "fileNameHandling.jl" begin

    # Need to read through these to see if we think this is 
    # the optimal caching strategy

    rtest = (10, "dr17", "dr17", "lco25m", "000+04", "10155", "58656", 150)
    itest = (10, "dr17", "dr17", "apo25m0000010i", "180+60", "9667", "58126", 295)
    dtest = (10, "sdsswork/mwm", "daily", "apo25m", "101689", "6455", "59817", 295)

    @test cache_skyname(rtest[4:7]...) == "../local_cache/58656/sky_lco25m_000+04_10155_58656.jdat"
    @test cache_skyname(itest[4:7]...) == "../local_cache/58126/sky_apo25m0000010i_180+60_9667_58126.jdat"
    @test cache_skyname(dtest[4:7]...) == "../local_cache/59817/sky_apo25m_101689_6455_59817.jdat"

    @test cache_skynameSpec(rtest[4:end]...) == "../local_cache/58656/450/skySpec_lco25m_000+04_10155_58656_150.jdat"
    @test cache_skynameSpec(itest[4:end]...) == "../local_cache/58126/295/skySpec_apo25m0000010i_180+60_9667_58126_295.jdat"
    @test cache_skynameSpec(dtest[4:end]...)== "../local_cache/59817/295/skySpec_apo25m_101689_6455_59817_295.jdat"

    @test cache_starname(rtest[4:end]...) == "../local_cache/58656/450/star_lco25m_000+04_10155_58656_150.jdat"
    @test cache_starname(itest[4:end]...) == "../inject_local_cache/58126/295/star_apo25m0000010i_180+60_9667_58126_295.jdat"
    @test cache_starname(dtest[4:end]...) == "../local_cache/59817/295/star_apo25m_101689_6455_59817_295.jdat"

    @test getUtahBase("dr17", "dr17") == "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/"
    @test getUtahBase("dr18", "dr18") == "/uufs/chpc.utah.edu/common/home/sdss/dr18/spectro/apogee/redux/dr18/"
    @test getUtahBase("sdsswork/mwm", "daily") == "/uufs/chpc.utah.edu/common/home/sdss/sdsswork/mwm/apogee/spectro/redux/daily/"

    @test build_platepath(rtest[2:7]...,"a") == "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit/lco25m/000+04/10155/58656/asPlate-a-10155-58656.fits"
    @test build_platepath(itest[2:7]...,"a") == "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit/apo25m/180+60/9667/58126/apPlate-a-9667-58126.fits"
    @test build_platepath(dtest[2:7]...,"a") == "/uufs/chpc.utah.edu/common/home/sdss/sdsswork/mwm/apogee/spectro/redux/daily/visit/apo25m/101689/6455/59817/apPlate-a-6455-59817.fits"

    tstfname = Downloads.download(sdss_public_tst*"asPlate-a-10155-58656.fits")
    @test getFramesFromPlate(tstfname) == ["30940037", "30940038", "30940039", "30940040", "30940041", "30940042", "30940043", "30940044"]

    @test build_framepath(rtest[2:4]...,rtest[7],"30940037","a") == "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/exposures/apogee-s/58656/as1D-a-30940037.fits"
    #these do not exist, but just testing the n/s switch
    @test build_framepath(itest[2:4]...,itest[7],"30940037","a") =="/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/exposures/apogee-n/58126/ap1D-a-30940037.fits"
    @test build_framepath(dtest[2:4]...,dtest[7],"30940037","a") == "/uufs/chpc.utah.edu/common/home/sdss/sdsswork/mwm/apogee/spectro/redux/daily/exposures/apogee-n/59817/ap1D-a-30940037.fits"

    ## visit and visit2cframe test
    vrname = build_visitpath(rtest[2:end]...)
    @test vrname == "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit/lco25m/000+04/10155/58656/asVisit-dr17-10155-58656-151.fits"
    viname = build_visitpath(itest[2:end]...)
    @test viname == "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit/apo25m/180+60/9667/58126/apVisit-dr17-9667-58126-006.fits"
    vdname = build_visitpath(dtest[2:end]...)
    @test vdname == "/uufs/chpc.utah.edu/common/home/sdss/sdsswork/mwm/apogee/spectro/redux/daily/visit/apo25m/101689/6455/59817/apVisit-daily-6455-59817-006.fits"

    @test visit2cframe(vrname,rtest[4],"30940037","a") == "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit/lco25m/000+04/10155/58656/asCframe-a-30940037.fits"
    @test visit2cframe(viname,itest[4],"30940037","a") == "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit/apo25m/180+60/9667/58126/apCframe-a-30940037.fits"
    @test visit2cframe(vdname,dtest[4],"30940037","a") == "/uufs/chpc.utah.edu/common/home/sdss/sdsswork/mwm/apogee/spectro/redux/daily/visit/apo25m/101689/6455/59817/apCframe-a-30940037.fits"
end