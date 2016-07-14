
# interatively get input parameters
while True:
    # WFDEI path:
    while True:
        input_path = raw_input("Please input the path for WFDEI data: ")
        input_path = os.path.realpath(os.path.expanduser(input_path))
        print input_path
        if os.path.lexists(input_path):
            break
        else:
            print "No such directory. Try again..."

    # output path:
    while True:
        output_path = raw_input("Please input the path for output: ")
        output_path = os.path.realpath(os.path.expanduser(output_path))
        print output_path
        if os.path.lexists(output_path):
            break
        else:
            print "No such directory. Try again..."

    # year range:
    while True:
        year_start = int(raw_input(
            "Please input the start year (YYYY): "))
        year_end = int(raw_input(
            "Please input the end year (YYYY): "))
        print(1979 <= year_start <= year_end <= 2014)
        if 1979 <= year_start <= year_end <= 2014:
            break
        else:
            print "Please input valid years. Try again..."

    # coordinates:
    while True:
        lat = float(raw_input(
            "Please input the latitude (in deg): "))
        lon = float(raw_input(
            "Please input the longitude (in deg): "))
        print(-90 < lat < 90 and -180 < lon < 180)
        if -90 < lat < 90 and -180 < lon < 180:
            break
        else:
            print "Please input valid coordinates. Try again..."

    start = time.time()
    write_SUEWS_forcing_1h(input_path, output_path,
                           year_start, year_end, lat, lon)
    end = time.time()
    print('time used in processing:' + '%.2f' % (end - start) + ' s')

    t = raw_input('Do you want to quit? Y/N')
    if t == 'Y' or t == 'y':
        ftp.quit()
        break
