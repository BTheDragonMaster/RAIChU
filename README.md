# RAIChU: Reaction Analysis through Imaging of Chemical Units
## Python-based informatics  tool for the visualisation of natural product chemistry

RAIChU is a chemoinformatics tool designed to illuminate the intricate processes of scaffold diversification in natural product chemistry through the implementation of 34 tailoring reactions. These reactions encompass family-specific posttranslational modifications across 25 RiPP families, enabling the exploration of advanced RiPP chemistry. A unique feature of RAIChU is its capability to detect and visually highlight potential sites for these tailoring reactions on molecules, a critical step due to the current limitations in predicting the regio- and stereochemistry of tailoring enzymes.

Our tool is particularly adept at visualizing the sequence of tailoring reactions in both modular systems, where reactions occur post-NRP/PK core release, and in RiPP and terpene systems, where tailoring, cyclisation, and proteolytic cleavage can occur in any order, albeit sequentially within each category. RAIChU allows for the flexible ordering of reactions, accommodating various biosynthetic strategies, and can automatically generate visual biosynthetic models showcasing the complexity of natural product biosynthesis.

However, it's important to note that while RAIChU strives for accuracy, the nature of some reactions, such as cationic cascade reactions in terpenes, necessitates programmatic shortcuts, potentially leading to simplified representations of reaction intermediates. To address this, RAIChU offers a summarized pathway view, highlighting only the initial and final intermediates for a clearer, albeit simplified, overview.

View the general.py script to learn about RAIChU's functionalities

Explore our wiki to dive deeper into RAIChU's features and capabilities:
[Home](https://github.com/SophieVromans/RAIChU/wiki)

[Getting Started with RAIChU](https://github.com/SophieVromans/RAIChU/wiki/Getting-Started)

[Visualizing Tailoring Reactions and instructions on tailoring enzymes](https://github.com/SophieVromans/RAIChU/wiki/Tailoring-enzymes)

[Examples of implemented RiPP families](https://github.com/SophieVromans/RAIChU/wiki/Examples-of-RiPP-families-that-can-be-implemented)
