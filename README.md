# GaSP: Train and Apply a Gaussian Stochastic Process Model

GaSP R package, created by [William J. Welch](https://www.stat.ubc.ca/~will/) and [Yilin Yang](https://www.linkedin.com/in/yilinyang-/).  
See the documentation for the basic outline and some simple examples of GaSP functions, and see the vignette for a more detailed description of GaSP as well as some noteworthy implementation choices made by the authors.

## Changelogs:
* Version 1.0.1: 
    * Added a vignette for GaSP.
    * Fixed memory leak in C functions.
    * Minor bug fixes for error matrix console output and Fit C initialization when 'random_error = TRUE'.
    
* Version 1.0.2: 
    * PROTECT bugs fixed
    * Compilation warnings about function prototypes, declarations, arguments fixed

* Version 1.0.3: 
    * More C compilation warnings fixed
    * R class() comparison with string fixed 

* Version 1.0.4
    * sprintf and vsprintf replaced by snprintf and vsnprintf, respectively

* Version 1.0.5
    * C format specifiers and type casts fixed for output messages
    
* Version 1.0.6
    * C types and type casts fixed