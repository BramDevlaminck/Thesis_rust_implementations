### Requirements
`clang` and `cmake` should be installed on the machine since it is used under the hood to create the C bindings in [`libsais64-rs`](libsais64-rs).  
On Debian Linux this can be installed by executing
```sh
 sudo apt install clang
```
```sh
 sudo apt install cmake
```
On MacOS `clang` is part of the Xcode Command Line Tools.
`cmake` can be installed using [brew](https://brew.sh/)
```sh
brew install cmake
```

### Installation
Clone this repository with the following command:
```sh
https://github.com/BramDevlaminck/Thesis_rust_implementations.git
```
and initialize the git submodules by executing
```sh
git submodule update --init --recursive
```