import os
import sys


# TODO: Duplicated in Gary's
def tpad(time, length=7):

    # TODO: Check that time is valid
    time = list(time)

    assert(len(time) > 2)

    if len(time) > length:
        time = time[0:length]
    else:
        pad = length - len(time)
        time = time + pad*[0]

    return tuple(time)


# TODO: Duplicated in Gary's
def tstr(time, length=7):
    """Create date/time string of the convention to tag files with given array of integers

    tstr((2000, 1, 1, 2)) # 2000:01:01T02:00:00
    tstr((2000, 1, 1, 2, 3)) # 2000:01:01T02:03:00
    tstr((2000, 1, 1, 2, 3, 4)) # 2000:01:01T02:03:04
    tstr((2000, 1, 1, 2, 3, 4, 567)) # 2000:01:01T02:03:04.567
    """
    import datetime

    if isinstance(time, datetime.date ):
        return time.strftime("%Y-%m-%dT%H:%M:%SZ")

    else:
        # ISO 8601
        assert(len(time) > 2)

        if length == 7:
            return '%04d-%02d-%02dT%02d:%02d:%02d.%03d' % tpad(time, length=length)
        elif length == 6:
            return '%04d-%02d-%02dT%02d:%02d:%02d' % tpad(time, length=length)
        elif length == 5:
            return '%04d-%02d-%02dT%02d:%02d' % tpad(time, length=length)
        elif length == 4:
            return '%04d-%02d-%02dT%02d' % tpad(time, length=length)
        elif length == 3:
            return '%04d-%02d-%02d' % tpad(time, length=length)

# TODO: Duplicated in Gary's
def time2datetime(t):
    import datetime as dt

    for i in range(len(t)):
        if int(t[i]) != t[i]:
            raise ValueError("int(t[{0:d}] != t[{0:d}] = {1:f}".format(i, t[i]))\

    if len(t) < 3:
        raise ValueError('Time list/tuple must have 3 or more elements')
    if len(t) == 3:
        return dt.datetime(int(t[0]), int(t[1]), int(t[2]))
    if len(t) == 4:
        return dt.datetime(int(t[0]), int(t[1]), int(t[2]), int(t[3]))
    if len(t) == 5:
        return dt.datetime(int(t[0]), int(t[1]), int(t[2]), int(t[3]), int(t[4]))
    if len(t) == 6:
        return dt.datetime(int(t[0]), int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5]))
    if len(t) == 7:
        return dt.datetime(int(t[0]), int(t[1]), int(t[2]), int(t[3]), int(t[4]), int(t[5]), int(t[6]))


def prompt(question, default=''):
    # Based on suggestions in https://gist.github.com/garrettdreyfus/8153571

    import sys
    if sys.version_info[0] > 2:
        reply = str(input(question + ': ')).lower().strip()
    else:
        reply = str(raw_input(question + ': ')).lower().strip()

    if len(reply) == 0:
        return default
    return reply[:1]


def install_paraview(paraview_version, install_path='/tmp/'):


    def extract(fullpath, install_path):
        cmd_list = ['tar', 'zxf', '--checkpoint', 10000, '--checkpoint-action="echo=Extracted file %u/124000"',fullpath, '--directory', install_path]
        print('Executing ' + ' '.join(cmd_list))
        import subprocess
        subprocess.run(cmd_list)

    file = paraview_version + '.tar.gz'
    tmpdir = '/tmp/'
    fullpath = tmpdir + file

    if os.path.exists(fullpath):
        extract(fullpath, install_path)
        return

    print('Downloading ' + file + ' to ' + tmpdir)
    url = 'https://www.paraview.org/paraview-downloads/download.php?' + \
          'submit=Download&version=v5.8&type=binary&os=Linux&downloadFile=' + file

    # https://stackoverflow.com/questions/15644964/python-progress-bar-and-downloads
    from clint.textui import progress
    import requests

    r = requests.get(url, stream=True)
    with open(fullpath, 'wb') as f:
        total_length = int(r.headers.get('content-length'))
        for chunk in progress.bar(r.iter_content(chunk_size=1024), expected_size=(total_length/1024) + 1):
            if chunk:
                f.write(chunk)
                f.flush()


def compatability_check(debug=False):

    import os
    import glob
    import site
    import warnings
    import subprocess

    debug = True
    if debug:
        print('sys.platform = ' + sys.platform)
    if sys.platform.startswith("darwin"):
        versions = glob.glob("/Applications/ParaView*")
        versions.sort()
        version = versions[-1]
        if len(version) == 0:
            print("ParaView not found in /Applications directory." \
                " To install ParaView, see https://www.paraview.org/download/.")
            # TODO?: Allow for ParaView not located in /Applications/
            sys.exit(1)
        use = ""
        # versions list will have elements of, e.g.,
        # ["/Applications/ParaView-5.7.0.app", "/Applications/ParaView-5.8.0.app"].
        for version in versions:
            if debug:
                print("Found " + version)
            version_str = version \
                            .replace("/Applications/ParaView-", "") \
                            .replace(".app", "")
            version_num = [int(i) for i in version_str.split(".")]
            # version_num will be tuple of (major, minor, patch)
            if list(sys.version_info[0:2]) > [2, 7] and version_num <= [5, 8]:
                if debug:
                    print('Python > 2.7 is not compatible with ParaView ' \
                        + version_str + ' on OS-X.')
            elif list(sys.version_info[0:2]) == [2, 7] and version_num[0:2] == [5, 8]:
                if debug:
                    print('A bug in ParaView 5.8 for OS-X prevents magnetovis from working. ' \
                        + 'See https://gitlab.kitware.com/paraview/paraview/-/issues/20146')
            else:
                use = version
                PARAVIEW = use + "/Contents/MacOS/paraview"
                PVPYTHON = use + "/Contents/bin/pvpython"
        if use == "":
            sys.exit(1)
        else:
            if debug:
                print("Using " + use)

    elif sys.platform.startswith("linux"):

        if sys.version_info[0] < 3:
            from backports import configparser
        else:
            import configparser

        config_file = os.path.expanduser('~') + '/.magnetovis.conf'
        PARAVIEWPATH = ''
        if 'PARAVIEWPATH' in os.environ:
            if debug:
                print('PARAVIEW os environment variable set as ' + os.environ['PARAVIEWPATH'])
            PARAVIEWPATH = os.path.expanduser(config['DEFAULT']['PARAVIEWPATH'])
            PARAVIEW = PARAVIEWPATH + '/bin/paraview'
            PVPYTHON = PARAVIEWPATH + '/bin/pvpython'
        else:
            if debug:
                print('PARAVIEW os environment not set.')

        if PARAVIEWPATH == '' and os.path.exists(config_file):
            if debug:
                print('File ' + config_file + ' found')
            config = configparser.ConfigParser()
            config.read(config_file)
            if 'PARAVIEWPATH' in config['DEFAULT']:
                PARAVIEWPATH = os.path.expanduser(config['DEFAULT']['PARAVIEWPATH'])
                if debug:
                    print('Found PARAVIEWPATH = ' + PARAVIEWPATH + ' in [DEFAULT] section of ' + config_file + '.')
                if os.path.exists(PARAVIEWPATH):
                    PARAVIEW = PARAVIEWPATH + '/bin/paraview'
                    PVPYTHON = PARAVIEWPATH + '/bin/pvpython'
                else:
                    print(PARAVIEWPATH + ' not found.')
                    PARAVIEWPATH = ''
                    #sys.exit(1)
            else:
                if debug:
                    print('PARAVIEWPATH variable not found in ' + config_file + '.')

        if PARAVIEWPATH == '':
            PVPYTHON = ''
            try:
                if debug:
                    print('Trying output of `which paraview`')
                PARAVIEW = subprocess.check_output(['which','paraview'])
                try:
                    if debug:
                        print('Trying output of `which pvpython`')
                    PVPYTHON = subprocess.check_output(['which','pvpython'])
                except:
                    print("Executable named 'pvpython' was not found in path, but 'paraview' found at " \
                          + PARAVIEW + ". Check installation; pvpython binary should be in the " \
                          + "same directory as paraview binary.")
                    sys.exit(1)
            except:
                print("Executable named 'paraview' was not found in path.")
                if prompt('Attempt to install? [y]/n', default='y') == 'y':
                    install_path = prompt('Enter installation path [~]', default='~')
                    install_path = os.path.expanduser(install_path)
                    paraview_version = 'ParaView-5.8.0-MPI-Linux-Python3.7-64bit'
                    install_paraview(paraview_version, install_path=install_path)
                    PARAVIEWPATH = install_path + '/' + paraview_version
                    PARAVIEW = PARAVIEWPATH + '/bin/paraview'
                    PVPYTHON = PARAVIEWPATH + '/bin/pvpython'
                    config = configparser.ConfigParser()
                    if os.path.exists(config_file):
                        if debug:
                            print("Found " + config_file + ". Updating PARAVIEWPATH " + \
                                  "variable to " + PARAVIEWPATH)
                        config.read(config_file)
                        if 'PARAVIEWPATH' in config['DEFAULT']:
                            config['DEFAULT']['PARAVIEWPATH'] = PARAVIEWPATH
                        else:
                            config['DEFAULT'] = {'PARAVIEWPATH': PARAVIEWPATH}
                        with open(config_file, 'w') as f:
                            config.write(f)
                        if debug:
                            print("Updated PARAVIEWPATH variable to " + PARAVIEWPATH)
                    else:
                        if debug:
                            print("Writing " + config_file)
                        config['DEFAULT'] = {'PARAVIEWPATH': install_path}
                        with open(config_file, 'w') as f:
                            config.write(f)
                        if debug:
                            print("Wrote " + config_file)
                else:
                    sys.exit(1)

    else:
        print("Installation implemented only for Linux and OS-X.")
        sys.exit(0)

    util_path = os.path.dirname(os.path.realpath(__file__))
    ver_path = os.path.join(util_path,'..','etc','version.py')
    if debug:
        print("Executing %s %s" % (PVPYTHON, ver_path))
    try:
            pvpython_version_str = subprocess.check_output([PVPYTHON, ver_path])
    except:
        print("Could not execute " + PVPYTHON + " " + ver_path + ". Exiting.")
        sys.exit(1)

    ver = sys.version_info
    python_version_str = str(ver.major) + "." + str(ver.minor) + "." + str(ver.micro)

    pvpython_version_str = pvpython_version_str.decode().strip()
    pvpython_version_lst = [int(i) for i in pvpython_version_str.split(".")]

    if pvpython_version_lst[0:2] == [2, 7] and sys.version_info[0:2] != (2, 7):
        print('Installed ParaView uses Python 2.7 and so magnetovis requires Python 2.7. Exiting.')
        print("pvpython Python version           = " + pvpython_version_str)
        print("Python version executing setup.py = " + python_version_str)
        sys.exit(1)

    if sys.version_info[0:2] != pvpython_version_lst[0:2]:
        # Hack to prevent warning from including unneeded information.
        def _warning(message, category=UserWarning, filename='', lineno=-1, file=None, line=''):
            print("UserWarning: " + str(message))
        warnings.showwarning = _warning

        warnings.warn("Installed ParaView uses Python " + pvpython_version_str \
                    + " and magnetovis is being installed for Python " \
                    + python_version_str \
                    + ". There may be compatability issues if you use " \
                    + "libraries and code used are not compatable with Python " \
                    + pvpython_version_str + ".")

    return PARAVIEW
