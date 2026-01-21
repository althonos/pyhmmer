if TARGET_SYSTEM == "Linux":
    from .cookie cimport (
        _LinuxReader as _FileobjReader,
        _LinuxWriter as _FileobjWriter,
    )
elif TARGET_SYSTEM == "Windows":
    from .win32 cimport (
        _Win32SynchronizedReader as _FileobjReader, 
        _Win32SynchronizedWriter as _FileobjWriter
    )
elif TARGET_SYSTEM == "Darwin" or TARGET_SYSTEM.endswith("BSD"):
    from .fun cimport (
        _BSDReader as _FileobjReader, 
        _BSDWriter as _FileobjWriter
    )
