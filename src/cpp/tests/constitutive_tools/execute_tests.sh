make && ./test_constitutive_tools && cat results.tex && cat results.tex | grep -i False || printf "\nNo failures\n"
