if TARGET_SYSTEM == "Linux":
    from .cookie cimport (
        fileobj_linux_open as fopen_obj,
        _LinuxReader as _FileobjReader,
        _LinuxWriter as _FileobjWriter,
    )
elif TARGET_SYSTEM == "Windows":
    from .win32 cimport (
        fopen_obj,
        _Win32SynchronizedReader as _FileobjReader, 
        _Win32SynchronizedWriter as _FileobjWriter
    )
elif TARGET_SYSTEM == "Darwin" or TARGET_SYSTEM.endswith("BSD"):
    from .fun cimport (
        fileobj_bsd_open as fopen_obj,
        _BSDReader as _FileobjReader, 
        _BSDWriter as _FileobjWriter
    )
