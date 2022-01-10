# MGP code repo

In this repository, you will find the code related to the eLife publication "[Relating multivariate shapes to genescapes using phenotype-biological process associations for craniofacial shape](https://elifesciences.org/articles/68623)".

The root directory is structured primarily as a Shiny app, although you can find the scripts used for analysis (and revisions) under the [analyses](/analyses) folder. In that folder you will also find an example script for using the API from within R. Finally, the code that defines the API can be found under the [MGP_API](/MGP_API) folder. If you're curious about how the app/API are deployed at [genopheno.ucalgary.ca/MGP](genopheno.ucalgary.ca/MGP), have a look at the various [configuration](/configurations) files for the shiny server, plumber service, and nginx reverse proxy.

If you want to learn Shiny, come to my workshop with [CalgaryR](https://www.meetup.com/calgaryr/) meetup group on December 8th, 2021. If you missed that date, the repo is linked [here](https://github.com/J0vid/CalgaryR_Shiny).
