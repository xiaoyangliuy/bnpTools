# Documentation to use bnpTools

## Using APSshare installed packages and adding bnpTools in python path for execution

> Using installed python packages in APSshare: 
>  - Open terminal 
    - Type "bash"
    - If the prefer coding enviornment is Jupyter notebook, below is the command to use: "/APSshare/anaconda3/x86_64/bin/jupyter notebook"
    - Usuually it will either append a new browser tab or open up a new FireFox session for jupyter notebook 
    - If none of the above happen, we can manually open a new page of FireFox, and type "localhost:xxxx" in the url tab, where 'xxxx' is the port number, which can be found in Terminal output after jupyter notebook command is executed
    
>    - To enable remote coding using Jupyter notebook, we can use ssh tunneling to transmit data by specifying a listening port upon starting the jupyter notebook
    > - ssh to the computer with defined listening port: <br>
    ssh -t -t userbnp@ivy.xray.aps.anl.gov -L 1980:localhost:1980
    - Execute jupyter notebook with specified listening port: <br>
    /APSshare/anaconda3/x86_64/bin/jupyter notebook --no-browser --port=1987 --allow-root --ip=0.0.0.0
    - In the local computer, copy and paste the url that links to the notebook on the remote machine (ivy). The url usually starts with "127.0.0. ...."  
    
> Adding bnpScan in PATH: <br>


```python
import sys
print(sys.path)
# sys.path.append('/home/beams/USERBNP/scripts/bnpTools')
# from bnpScan import bnpScan
```

    ['/home/beams11/USERBNP/scripts/bnpTools', '/APSshare/anaconda3/x86_64/lib/python37.zip', '/APSshare/anaconda3/x86_64/lib/python3.7', '/APSshare/anaconda3/x86_64/lib/python3.7/lib-dynload', '', '/APSshare/anaconda3/x86_64/lib/python3.7/site-packages', '/APSshare/anaconda3/x86_64/lib/python3.7/site-packages/locket-0.2.1-py3.7.egg', '/APSshare/anaconda3/x86_64/lib/python3.7/site-packages/pyRestTable-2020.0.2-py3.7.egg', '/APSshare/anaconda3/x86_64/lib/python3.7/site-packages/IPython/extensions', '/home/beams11/USERBNP/.ipython', '/home/beams/USERBNP/scripts/bnpTools']



```python

```
