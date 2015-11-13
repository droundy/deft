# -*- mode: ruby -*-
# vi: set ft=ruby :

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
  # If true, then any SSH connections made will enable agent forwarding.
  # Default value: false
  # config.ssh.forward_agent = true

  config.vm.define "trusty" do |trusty|
    trusty.vm.box = "ubuntu/trusty"

    trusty.vm.provider "virtualbox" do |v|
      v.memory = 16*1024
      v.cpus = 2
    end
    trusty.vm.provision :shell, path: "vagrant/provision-debian.sh"
  end

  config.vm.define "jessie" do |jessie|
    jessie.vm.box = "debian/jessie64"

    jessie.vm.provider "virtualbox" do |v|
      v.memory = 16*1024
      v.cpus = 2
    end
    jessie.vm.provision :shell, path: "vagrant/provision-debian.sh"
  end

  # config.vm.define "freebsd64" do |freebsd64|
  #   freebsd64.vm.box = "chef/freebsd-10.0"
  #   #freebsd64.vm.autostart = false

  #   freebsd64.vm.provision :shell, path: "vagrant/provision-freebsd.sh"
  # end

  # config.vm.define "win7" do |win7|
  #   win7.vm.box = "http://aka.ms/vagrant-win7-ie11"
  #   win7.vm.autostart = false
  # end
end
