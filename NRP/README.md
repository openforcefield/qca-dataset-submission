## When running a job

Do not waste resources. Our access to NRP is not guaranteed and is dependent on our good behavior. If we get in trouble they will apply penalties to the entire openforcefield namespace.

* Have a browser tab open to monitor your % utilization. Keep this as high as possible. Generally, GPU+processor requests should be over 70% and memory over 35%
* Have the element/matrix chat open and check every few hours whether the admins are trying to get your attention.
* See the QCF manager yaml for some of the special keywords/tolerations that make workers pre-emptible. 

# Playbook

Tolerations are special resources that are preemptible

If changing manager/deployment CPUs or RAM, coordinate changes in both manager.yaml and deployment.yaml, and then delete and remake appropriate manager.yaml secret

Deployments have tons of tolerations allowing them to run on pre-emptible resources that normal jobs can’t use, also have nice priority. This means our jobs will die more often, but we should be comfortable requesting tons of resources. Just make sure we use what we request, or scale down if not!

    kubectl --context=nautilus create secret generic manager-qm-config-yaml --from-file=manager.yaml

    kubectl --context=nautilus apply -f deployment.yaml

    kubectl --context=nautilus get pods

    kubectl logs openff-qca-qm-856ccd6795-pnkkc

    kubectl --context=nautilus delete deployment openff-qca-qm

    kubectl --context=nautilus delete secret manager-qm-config-yaml

    kubectl --context=nautilus get deployment openff-qca-qm

    kubectl scale --replicas=25 -f deployment.yaml

# Monitoring

https://grafana.nrp-nautilus.io/dashboards

In particular 

    Utilization: https://grafana.nrp-nautilus.io/d/zytktgwWz/utilization-by-namespace?orgId=1&from=now-3h&to=now
    Bacon plots: https://grafana.nrp-nautilus.io/d/85a562078cdf77779eaa1add43ccec1e/kubernetes-compute-resources-namespace-pods?orgId=1&refresh=10s&var-datasource=default&var-cluster=&var-namespace=openforcefield&from=now-30m&to=now

# Onboarding/getting started

### Making account on NRP

* Make an NRP account using your @openforcefield.org email. https://docs.nationalresearchplatform.org/userdocs/start/get-access/
* After making an NRP account, sign into matrix chat by making (yet another) new account. It's important to monitor the chat while running jobs.
* Have the openforcefield admin add you to the namespace. If Jeff has been hit by a bus ask on the support chat to have someone new promoted to admin.
* Click "Get Config" on the top right of the login page and follow the instructions.
* Install kubectl on your computer (brew using mac, several options on linux) https://kubernetes.io/docs/tasks/tools/#kubectl
* Verify that your config is happy by running something like `kubectl --context=nautilus get deployment` and ensure you see OpenFF's running jobs. 

### Getting QCA credentials

JW will share worker credentials via LastPass. 

### Your first job

* Copy `qca-dataset-submission/NRP/manager_template.yaml` and `deployment_template.yaml` to `manager.yaml` and `deployment.yaml` in a new folder for your personal use.
* Modify `manager.yaml`:
    * Substitute the QCA username and password into the manager.yaml file.
    * Change the cores_per_worker or memory_per_worker in manager.yaml to something between 1 and 10, and 6 and 30 respectively
* Create a new secret with your `manager.yaml`, appending your initials to the secret name (eg, mine would be `manager-qm-config-yaml-jw`)
* Update `deployment.yaml`:
    * append your initials to all instances of `openff-qca-qm` (eg, mine would be `openff-qca-qm-jw`)
    * update the cpu and memory requests at the bottom to be equal to the cores and memory you set in `manager.yaml` (sometimes we have to add some buffer for memory since psi4/qcengine exceeds the limits in `manager.yaml`)
* Launch your jobs using `kubectl apply ...`
* Get your pod names using eg. `kubectl --context=nautilus get pods | grep openff-qca-qm-jw`
* Check out some logs using `kubectl logs <podname>`
* Find your jobs on [grafana](https://grafana.nrp-nautilus.io/d/85a562078cdf77779eaa1add43ccec1e/kubernetes-compute-resources-namespace-pods?orgId=1&refresh=10s&var-datasource=default&var-cluster=&var-namespace=openforcefield&from=now-30m&to=now)
* Increase your number of workers to 4 using `kubectl scale ...`
* Have fun, watch them go for a bit
* Shut down all your workers using `kubectl delete ...`

### Other topics

* To stay in NRP's good graces, we should update our publications and stuff in our namespace. Would anyone be up to give that a swing?
* This might be suitable for fitting runs? There's lots I don't know about having kubes talk to each other and persistent storage. And we'd want to remove tolerances/pre-emtibility from the deployment to ensure we don't lose info.

