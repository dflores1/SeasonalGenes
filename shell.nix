let
  sources = import ./nix/sources.nix;
  pkgs = import sources.nixpkgs {};
  # R and bioconductor pkgs are not fine-tuned for apple m1 processors. Use x86 in that case.
  change_to_x86 = with pkgs.stdenv.hostPlatform; isDarwin && isAarch64;
  R_pkgs = if change_to_x86 then (import sources.nixpkgs { system = "x86_64-darwin"; }) else pkgs;
  myR = with R_pkgs; rWrapper.override {
    packages = with rPackages; [
	ggplot2
	gridExtra
	gridtext
	ggpubr
	ggpp
	lubridate
	lme4
	pbapply
	stringr
	readxl
	HiClimR
	hgu133a_db
	reshape2
    ];
  };
in
  pkgs.mkShell {
    packages = [ pkgs.xz ]; # If you want an editor like emacs include here: pkgs.emacsPackages.emacs pkgs.emacsPackages.ess
    inputsFrom = [ myR ];
    shellHook = ''
    echo From here you can run your favourite editor. Or simply run the R script. Enjoy!
    '';
  }
