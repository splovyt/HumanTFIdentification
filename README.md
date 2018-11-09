# Li Lab UCSF - Human TF Identification

A web application that takes gene expression fold changes as input and returns human transcription factors crucial to the experimental condition.

The correlation score of the 'impact' of the TF regarding the changes in gene expression is returned together with the DNA-binding motif on which the score is based. The transcription factors are sorted by descending correlation score.

This application is made in Django, Python, HTML, CSS, and JavaScript (D3).

### Example output of an old-vs-young gene expression contrast:

![alt text](https://raw.githubusercontent.com/splovyt/HumanTFIdentification/master/web_interface_screenshot.png)
