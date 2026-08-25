{
    description = "An FHS shell for julia.";

    inputs = {
        nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
        # nixpkgs.url = "github:NixOS/nixpkgs/nixos-23.11";
    };

    outputs = { self, nixpkgs, home-manager }:
        let
            pkgs = nixpkgs.legacyPackages.x86_64-linux;
        in {
            devShell.x86_64-linux = (pkgs.buildFHSEnv {
                name = "conda";
                targetPkgs = pkgs: (
                    with pkgs; [
                        desktop-file-utils
                        libxcomposite
                        libxtst
                        libxrandr
                        libxext
                        libx11
                        libxfixes
                        libGL

                        gst_all_1.gstreamer
                        gst_all_1.gst-plugins-ugly
                        gst_all_1.gst-plugins-base
                        libdrm
                        xkeyboardconfig
                        libpciaccess

                        glib
                        gtk2
                        gtk3
                        bzip2
                        zlib
                        gdk-pixbuf

                        libxinerama
                        libxdamage
                        libxcursor
                        libxrender
                        libxScrnSaver
                        libxxf86vm
                        libxi
                        libSM
                        libICE
                        freetype
                        curlWithGnuTls
                        nspr
                        nss
                        fontconfig
                        cairo
                        pango
                        expat
                        dbus
                        cups
                        libcap
                        SDL2
                        libusb1
                        udev
                        dbus-glib
                        atk
                        at-spi2-atk
                        libudev0-shim

                        libxt
                        libxmu
                        libxcb
                        xcbutil
                        xcbutilwm
                        xcbutilimage
                        xcbutilkeysyms
                        xcbutilrenderutil
                        libGLU
                        libuuid
                        libogg
                        libvorbis
                        SDL
                        SDL2_image
                        glew110
                        openssl
                        libidn
                        tbb
                        wayland
                        xwayland
                        mesa
                        libxkbcommon
                        vulkan-loader

                        flac
                        freeglut
                        libjpeg
                        libpng12
                        libpulseaudio
                        libsamplerate
                        libmikmod
                        libtheora
                        libtiff
                        pixman
                        speex
                        SDL_image
                        SDL_ttf
                        SDL_mixer
                        SDL2_ttf
                        SDL2_mixer
                        libappindicator-gtk2
                        libcaca
                        libcanberra
                        libgcrypt
                        libvpx
                        librsvg
                        libxft
                        libvdpau
                        alsa-lib

                        harfbuzz
                        e2fsprogs
                        libgpg-error
                        keyutils.lib
                        libjack2
                        fribidi
                        p11-kit

                        gmp

                        # libraries not on the upstream include list, but nevertheless expected
                        # by at least one appimage
                        libtool.lib # for Synfigstudio
                        libxshmfence # for apple-music-electron
                        at-spi2-core
                        # julia
                        # libx11
                        # libxcb
                        # libxcomposite
                        # libxcursor
                        # libxdamage
                        # libxext
                        # libxfixes
                        # libxi
                        # libxrender
                        # libxtst
                        # libxrandr
                        # libxScrnSaver
                        # xwayland
                        # alsa-lib
                        # at-spi2-core
                        # dbus
                        # pango
                        # expat
                        # mesa
                        # cups
                        # libdrm
                        # gcc
                        # nspr
                        # nss
                        # cairo
                        # glib
                        # gtk3
                        # gdk-pixbuf
                    ]
                );
              profile = "";
            }).env;
        };
}
