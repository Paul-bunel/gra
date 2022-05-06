from urllib import request
from bs4 import BeautifulSoup
import re

headers = []

for i in range(1, 3068):
    # Fetch the html file
    url = 'https://bitterdb.agri.huji.ac.il/Receptor.php?id=' + str(i)
    # print(url)
    try:
        response = request.urlopen(url)
        html_doc = response.read()

        # Parse the html file
        soup = BeautifulSoup(html_doc, 'html.parser')

        for div in soup.find_all('div'):
            if 'id' in div.attrs and div.attrs['id'] == "curRecepSequence":
                prot_seq = div

        prot_seq = prot_seq.string
        prot_seq = re.sub('(\t|\r| +)', '', prot_seq)

        title = soup.find('h2')
        header = re.sub('(\n|\t|\r+)', '', title.contents[0])

        fasta_prot = '>' + header + prot_seq + '\n'

        # print(repr(fasta_prot))

        with open("test.fasta", 'a') as file:
            if header not in headers:
                file.write(fasta_prot)
        headers.append(header)
    except:
        # print(i)
        continue
