### Text diagrams extraction

Such feature only available from MadGraph5 versions higher than 3.8.0 (see pr [#398](https://github.com/mg5amcnlo/mg5amcnlo/pull/398) for more details)

After generating and outputting the process, use:
```bash
display diagrams_text ./diagrams/proc_name --no_open --merge
```
to extract and merge all the text diagrams. You might want give the merged file a different name other than `all_diagrams_text.txt`. You can do this inside the MadGraph5 CLI via
```bash
shell mv ./diagrams/proc_name/all_diagrams_text.txt ./diagrams/proc_name/all_diagrams_proc_name.txt
```
or outside the MG5 CLI using same command without the `shell`.
