{
  description = "Development environment for Fortran with BLAS";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  };
  outputs = { self, nixpkgs }: 
    let
      system = "x86_64-linux";
      
      pkgs = import nixpkgs { 
        inherit system;
      };
    in
    {
    devShells.x86_64-linux.default = pkgs.mkShell {
      buildInputs = [ pkgs.gfortran9 pkgs.blas ];
      shellHook = 
        ''
        '';
      BUILD_COMMAND = "make build";
      TEST_COMMAND  = "make run";
      CLEAN_COMMAND = "make clean";
    };
  };
}
